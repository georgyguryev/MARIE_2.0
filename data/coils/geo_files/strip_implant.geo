
 P1101 = 1101;
 Point(P1101) = {0.002500000000000002, -0.002500000000000002, 0.002500000000000002, 0.005000000000000004};
 P1102 = 1102;
 Point(P1102) = {0.002500000000000002, 0.002500000000000002, 0.002500000000000002, 0.005000000000000004};
 P1103 = 1103;
 Point(P1103) = {0.002500000000000002, -0.002500000000000002, -0.04750000000000004, 0.005000000000000004};
 P1104 = 1104;
 Point(P1104) = {0.002500000000000002, 0.002500000000000002, -0.04750000000000004, 0.005000000000000004};
 P1105 = 1105;
 Point(P1105) = {-0.007500000000000007, -0.002500000000000002, 0.002500000000000002, 0.005000000000000004};
 P1106 = 1106;
 Point(P1106) = {-0.007500000000000007, 0.002500000000000002, 0.002500000000000002, 0.005000000000000004};
 P1107 = 1107;
 Point(P1107) = {-0.007500000000000007, -0.002500000000000002, 0.03750000000000003, 0.005000000000000004};
 P1108 = 1108;
 Point(P1108) = {-0.007500000000000007, 0.002500000000000002, 0.03750000000000003, 0.005000000000000004};

  

  
 L1501 = 1501;
 Line(L1501) = {P1103, P1101};
 L1502 = 1502;
 Line(L1502) = {P1101, P1102};
 L1503 = 1503;
 Line(L1503) = {P1102, P1104};
 L1504 = 1504;
 Line(L1504) = {P1104, P1103};
 L1505 = 1505;
 Line(L1505) = {P1101, P1105};
 L1506 = 1506;
 Line(L1506) = {P1105, P1106};
 L1507 = 1507;
 Line(L1507) = {P1106, P1102};
 L1508 = 1508;
 Line(L1508) = {P1105, P1107};
 L1509 = 1509;
 Line(L1509) = {P1107, P1108};
 L1510 = 1510;
 Line(L1510) = {P1108, P1106};
 O1701 = 1701;
 Line Loop(O1701) = { L1501, L1502, L1503, L1504 };
 O1702 = 1702;
 Line Loop(O1702) = { L1505, L1506, L1507, -L1502 };
 O1703 = 1703;
 Line Loop(O1703) = { L1508, L1509, L1510, -L1506 };

  

  
 S1801 = 1801;
 Ruled Surface(S1801) = {O1701};
 S1802 = 1802;
 Ruled Surface(S1802) = {O1702};
 S1803 = 1803;
 Ruled Surface(S1803) = {O1703};

  
 Coil1 = 1000;
 Physical Surface(Coil1) = { S1801, S1802, S1803 };

 // --------------------------------------------------------------- 

