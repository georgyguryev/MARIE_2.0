% compile mex function assembly_lop_coup.cpp

%% TODO in the following releases

% try
% 	if ispc
% 		mex -output ompQuadCoil2Scat COMPFLAGS="$COMPFLAGS /openmp" -v -O -largeArrayDims QuadCoil2Scat_ompmex.cpp Coupling.cpp
% 		mex -output ompQuadScat2Coil COMPFLAGS="$COMPFLAGS /openmp" -v -O -largeArrayDims QuadScat2Coil_ompmex.cpp Coupling.cpp
% 	else
% 		mex -output ompQuadCoil2Scat COMPFLAGS='$COMPFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' CXXOPTIMFLAGS='$CXXOPTIMFLAGS -fopenmp' -v -O -largeArrayDims QuadCoil2Scat_ompmex.cpp Coupling.cpp
% 		mex -output ompQuadScat2Coil COMPFLAGS='$COMPFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' CXXOPTIMFLAGS='$CXXOPTIMFLAGS -fopenmp' -v -O -largeArrayDims QuadScat2Coil_ompmex.cpp Coupling.cpp
% 	end
% catch me
%     % use the sequential version
%     fprintf(1, '\n\n\n  Warning: no openmp compatible compiler\n  Mex functions will be executed in single core\n\n');
%     pause(3);
% end
% mex -output mexQuadCoil2Scat -v -O -largeArrayDims QuadCoil2Scat_mex.cpp Coupling.cpp
% mex -output mexQuadScat2Coil -v -O -largeArrayDims QuadScat2Coil_mex.cpp Coupling.cpp

%% temporary mexing of the coupling scripts

mex -output mexAssembly_Lop_coup -v -O -largeArrayDims assembly_lop_coup.cpp