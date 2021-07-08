#include <cstdint>
#include <cstdbool>
#include <cstdlib>
#include <cmath>
#include <omp.h>
#include "mex.h"
#include "matrix.h"

#include "EMMC.h"
#include "Assembly/VIE/SSQ.h"
#include "Assembly/VIE/VVQ.h"

template <typename F1,typename F2=F1>
inline vec3<F2> mxArrayToVec3(mxArray const* const mx,bool s = 1) {
	F1 const* p = (F1 const*)mxGetData(mx);
	return vec3<F2>((F2)p[0],(F2)p[s],(F2)p[2*s]);
}

template< class T >
constexpr bool is_floating_point_v = std::is_floating_point<T>::value;
template <typename F>
void mexPrintNum(F x,char a='\0') {
	if(is_floating_point_v<F>)
		mexPrintf("%+.3e%c",double(x),a);
	else
		mexPrintf("%+d%c",(int)x,a);
}

template <typename F,bool pwl,VIE::OpType op>
inline void paramMexFunction(F const freq,vec3<F> const& dr,Size3 const& Ndom,mxArray* Ko) {
	VKernel<F,pwl,op>* const Kr = reinterpret_cast<VKernel<F,pwl,op>*>(mxGetData(Ko));
	VKernel<F,pwl,op>* const Ki = reinterpret_cast<VKernel<F,pwl,op>*>(mxGetImagData(Ko));
	F const k0 = EM_2PI_C*freq;
	constexpr size_t NqF = 2;
	constexpr size_t NqM = 5;
	constexpr size_t NqN = 12;
	Size3 const NevalM = {5,5,5};
	Size3 const NevalN = {2,2,2};
	if(Ndom.x > NevalM.x || Ndom.y > NevalM.y || Ndom.z > NevalM.z) {
		//VVQ<F,pwl,NqF,op>(Ndom,Ndom,dr,k0,Kr,Ki);
		VVQFunctor<F,pwl,NqF,op>::compute(Ndom,Ndom,dr,k0,Kr,Ki);
	}
	if(Ndom.x > NevalN.x || Ndom.y > NevalN.y || Ndom.z > NevalN.z) {
		//VVQ<F,pwl,NqM,op>(Ndom,NevalM,dr,k0,Kr,Ki);
		VVQFunctor<F,pwl,NqM,op>::compute(Ndom,min(Ndom,NevalM),dr,k0,Kr,Ki);
	}
	SSQ<F,pwl,NqN,op>(k0,Ndom,dr,Kr,Ki);
}

template <typename F,bool pwl>
inline void paramMexFunction(char const op,F const freq,vec3<F> const& dr,Size3 const& Ndom,mxArray* Ko) {
	switch(op) {
		case 'n':
		case 'N':
			paramMexFunction<F,pwl,VIE::N>(freq,dr,Ndom,Ko);
		break;
		case 'k':
		case 'K':
			paramMexFunction<F,pwl,VIE::K>(freq,dr,Ndom,Ko);
		break;
		default:
			mexErrMsgIdAndTxt("marie:vie:assembly:char_match","Unable to match mode to known modes; got %s, expected 'N' or 'K' (case-insensitive)",op);
	}
}

template <typename F>
inline void paramMexFunction(bool pwl,char const op,F const freq,vec3<F> const& dr,Size3 const& Ndom,mxArray* Ko) {
	if(pwl) {
		paramMexFunction<F,true>(op,freq,dr,Ndom,Ko);
	}
	else {
		paramMexFunction<F,false>(op,freq,dr,Ndom,Ko);
	}
}

void ensureStackSize() {
	static bool failure = false;
	static bool success = false;
	if(success) {
		// if you change it after we check,
		// then that's on you
		return;
	}
	if(failure) {
		mexErrMsgIdAndTxt("marie:omp:gomploading","A previous failure was detected. Make sure to start MATLAB with\n\tOMP_STACKSIZE and LD_PRELOAD\nset as explained in the README");
	}
	char* stacksize = getenv("OMP_STACKSIZE");
	char* ldpreload = getenv("LD_PRELOAD");
	if((failure = ldpreload == NULL)) {
		mexErrMsgIdAndTxt("marie:omp:ldpreload","LD_PRELOAD must be set when starting MATLAB, and must contain libgomp.so");
	}
	if((failure = stacksize == NULL)) {
		mexErrMsgIdAndTxt("marie:omp:stacksize","OMP_STACKSIZE must be set when starting MATLAB\nsetenv(...) within MATLAB will not do.");
	}
	char const numspan[] = "0123456789";
	char const sizespan[] = "kKmMgG";
	size_t len = strlen(stacksize);
	bool issize = strchr(sizespan,stacksize[len-1]) != NULL;
	if((failure = len > 12)) {
		mexErrMsgIdAndTxt("marie:omp:stacksize","OMP_STACKSIZE is suspiciously long; maximum length of 12 is supported; found %u.",len);
	}
	if((failure = !issize && strchr(numspan,stacksize[len-1]) == NULL)) {
		mexErrMsgIdAndTxt("marie:omp:stacksize","OMP_STACKSIZE must end with a number or a size modifier (k,m,g); found %c",stacksize[len-1]);
	}
	size_t dstart = strcspn(stacksize,numspan);
	size_t dsize = strspn(stacksize,numspan);
	if((failure = dstart > 0 || dsize+issize < len)) {
		mexErrMsgIdAndTxt("marie:omp:stacksize","OMP_STACKSIZE must consist of digits, followed by an optional modifier");
	}
	long int stacksizeint = strtol(stacksize,NULL,10);
	if(issize) {
		switch(stacksize[len-1]) {
			case 'g':
			case 'G':
				stacksizeint <<= 10;
			case 'm':
			case 'M':
				stacksizeint <<= 10;
			case 'k':
			case 'K':
				stacksizeint <<= 10;
		}
	}
	if((failure = stacksizeint < (1 << 22))) {
		mexErrMsgIdAndTxt("marie:omp:stacksize","OMP_STACKSIZE is too small; we recommend setting it to 4m or greater");
	}
	success = true;
}
#include <sys/resource.h>
void checkrlimit() {
	struct rlimit rlim;
	getrlimit(RLIMIT_STACK,&rlim);
	if(rlim.rlim_max == RLIM_INFINITY)
		mexPrintf("rlimit(soft=%d)\n",rlim.rlim_cur);
	else
		mexPrintf("rlimit(soft=%d,hard=%d)\n",rlim.rlim_cur,rlim.rlim_max);
}

void mexFunction(int nlhs, mxArray* plhs[],int nrhs,const mxArray* prhs[]) {
	if(nrhs < 4 || nrhs > 5) {
		mexErrMsgIdAndTxt("marie:vie:Nop:invalid_nargin","Expecting 4-5 inputs; got %d",nrhs);
	}
	if(nlhs > 1) {
		mexErrMsgIdAndTxt("marie:vie:Nop:invalid_nargout","Expecting <= 1 outputs; got %d",nlhs);
	}
	if(!mxIsScalar(prhs[0]) || !mxIsChar(prhs[0])) {
		mexErrMsgIdAndTxt("marie:vie:assembly:invalid_argin","Input %d is not a char scalar");
	}
	if(!mxIsScalar(prhs[1]) || !(mxIsDouble(prhs[1]) || mxIsSingle(prhs[1]))) {
		mexErrMsgIdAndTxt("marie:vie:assembly:invalid_argin","Input %d is not a double scalar",1);
	}
	if(!mxIsDouble(prhs[2]) || (mxGetNumberOfElements(prhs[2]) != 3 && mxGetNumberOfElements(prhs[2]) != 1)) {
		mexErrMsgIdAndTxt("marie:vie:assembly:invalid_argin","Input %d is not a double triplet",2);
	}
	if(!mxIsDouble(prhs[3]) || mxGetNumberOfElements(prhs[3]) != 3) {
		mexErrMsgIdAndTxt("marie:vie:assembly:invalid_argin","Input %d is not a double triplet",3);
	}
	if(nrhs == 5 && !mxIsLogicalScalar(prhs[4])) {
		mexErrMsgIdAndTxt("marie:vie:assembly:invalid_argin","Input %d is not a logical scalar",4);
	}
	//ensureStackSize();
	//checkrlimit();
	bool aniso = mxGetNumberOfElements(prhs[2]) == 3;
	bool pwl = nrhs < 5 || mxIsLogicalScalarTrue(prhs[4]);
	using uint = unsigned int;
	char* const op = mxArrayToString(prhs[0]);
	if(op == NULL) {
		mexErrMsgIdAndTxt("marie:vie:assembly:char_parse","Unable to parse first input as char array (either out-of-memory or invalid type)");
	}
	bool isN = *op == 'N' || *op == 'n';
	// get domain size and allocate output tensors
	Size3 Ndom = mxArrayToVec3<double,uint>(prhs[3]);
	size_t const dims[5] = {pwl ? VIE::Ni<1> : VIE::Ni<0>,isN ? VIE::Ng<VIE::N> : VIE::Ng<VIE::K>,(size_t)Ndom.x,(size_t)Ndom.y,(size_t)Ndom.z};
	plhs[0] = mxCreateNumericArray(5,dims,mxGetClassID(prhs[1]),mxCOMPLEX);
	// dispatch based on class of first input, freq
	if(mxIsSingle(prhs[1])) {
		float const freq = *((float const*)mxGetData(prhs[1]));
		vec3<float> dr = mxArrayToVec3<double,float>(prhs[2],aniso);
		paramMexFunction<float>(pwl,*op,freq,dr,Ndom,plhs[0]);
	}
	else {
		double const freq = *((double const*)mxGetData(prhs[1]));
		vec3<double> dr = mxArrayToVec3<double>(prhs[2],aniso);
		paramMexFunction<double>(pwl,*op,freq,dr,Ndom,plhs[0]);
	}
}
