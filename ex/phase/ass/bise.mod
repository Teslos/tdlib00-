<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<AssociatedSolution class=phase id=L
  equilibrated=0 debug=0>
    <SimpleSolution class=phase>
      <components> Bi Se </components>
      <Reference class=FuncTpx>
        <func_x> +x(Bi)* </func_x>
        <species> Bi
        </species>
        <func_x> +x(Se)* </func_x>
        <species> Se
        </species>
      </Reference> 
    </SimpleSolution> 
    <SimpleSolution class=phase>
      <components> Bi Se Bi2 Bi3 Bi4 Se2 Se3 Se4 Se5 Se6 Se7 Se8 BiSe </components>
      <Reference class=FuncTpx>
        <func_x> +x(Bi)* </func_x>
        <species> Bi
          <compound_Tp class=func_Tp>
            <IVT_Tp class=func_Tp>
              (
              <coef id=a unknown=0> 202568 </coef> +
              <coef id=b unknown=0> -48.3484 </coef> *T+
              <coef id=c unknown=0> -20.6738 </coef> *T*log(T)+
              <coef id=d unknown=0> -0.000198868 </coef> *T^2+
              <coef id=e unknown=0> 7.66273e-08 </coef> *T^3+
              <coef id=f unknown=0> -1.39641e-11 </coef> *T^4+
              <coef id=g unknown=0> -1285.96 </coef> /T)
            </IVT_Tp> 
            <limit_Tp> T = {298.15, 3000} p = {0, inf} </limit_Tp>
            <IVT_Tp class=func_Tp>
              (
              <coef> 261551 </coef> +
              <coef> -273.236 </coef> *T+
              <coef> 7.15263 </coef> *T*log(T)+
              <coef> -0.00567066 </coef> *T^2+
              <coef> 2.04373e-07 </coef> *T^3+
              <coef> -3.89225e-12 </coef> *T^4+
              <coef> -2.31076e+07 </coef> /T)
            </IVT_Tp> 
            <limit_Tp> T = {3000, 10000} p = {0, inf} </limit_Tp>
          </compound_Tp> 
        </species>
        <func_x> +x(Se)* </func_x>
        <species> Se
          <compound_Tp class=func_Tp>
            <IVT_Tp class=func_Tp>
              (
              <coef> 233135 </coef> +
              <coef> -73.3768 </coef> *T+
              <coef> -14.9411 </coef> *T*log(T)+
              <coef> -0.00761606 </coef> *T^2+
              <coef> 1.35404e-06 </coef> *T^3+
              <coef> -1.02282e-10 </coef> *T^4+
              <coef> -90656.3 </coef> /T)
            </IVT_Tp> 
            <limit_Tp> T = {298.15, 1500} p = {0, inf} </limit_Tp>
            <IVT_Tp class=func_Tp>
              (
              <coef> 235900 </coef> +
              <coef> -43.4885 </coef> *T+
              <coef> -20.2094 </coef> *T*log(T)+
              <coef> -0.000760228 </coef> *T^2+
              <coef> 3.59295e-08 </coef> *T^3+
              <coef> -7.44689e-13 </coef> *T^4+
              <coef> -2.03599e+06 </coef> /T)
            </IVT_Tp> 
            <limit_Tp> T = {1500, 10000} p = {0, inf} </limit_Tp>
          </compound_Tp> 
        </species>
        <func_x> +x(Bi2)* </func_x>
        <species> Bi2
          <compound_Tp class=func_Tp>
            <IVT_Tp class=func_Tp>
              (
              <coef> 207362 </coef> +
              <coef> 1.56625 </coef> *T+
              <coef> -41.2243 </coef> *T*log(T)+
              <coef> 0.00551565 </coef> *T^2+
              <coef> -1.73281e-06 </coef> *T^3+
              <coef> 1.17723e-10 </coef> *T^4+
              <coef> 76718.4 </coef> /T)
            </IVT_Tp> 
            <limit_Tp> T = {298.15, 1500} p = {0, inf} </limit_Tp>
            <IVT_Tp class=func_Tp>
              (
              <coef> 226104 </coef> +
              <coef> -217.738 </coef> *T+
              <coef> -9.30681 </coef> *T*log(T)+
              <coef> -0.0160165 </coef> *T^2+
              <coef> 1.26719e-06 </coef> *T^3+
              <coef> -4.48408e-11 </coef> *T^4+
              <coef> -1.07796e+06 </coef> /T)
            </IVT_Tp> 
            <limit_Tp> T = {1500, 6000} p = {0, inf} </limit_Tp>
          </compound_Tp> 
        </species>
        <func_x> +x(Bi3)* </func_x>
        <species> Bi3
          <compound_Tp class=func_Tp>
            <IVT_Tp class=func_Tp>
              (
              <coef> 271073 </coef> +
              <coef> 20.1435 </coef> *T+
              <coef> -58.184 </coef> *T*log(T)+
              <coef> -1.86769e-05 </coef> *T^2+
              <coef> 5.02083e-09 </coef> *T^3+
              <coef> -6.98738e-13 </coef> *T^4+
              <coef> 18521.1 </coef> /T)
            </IVT_Tp> 
            <limit_Tp> T = {298.15, 1500} p = {0, inf} </limit_Tp>
            <IVT_Tp class=func_Tp>
              (
              <coef> 266889 </coef> +
              <coef> 45.756 </coef> *T+
              <coef> -61.5985 </coef> *T*log(T)+
              <coef> 0.00128258 </coef> *T^2+
              <coef> -1.06436e-07 </coef> *T^3+
              <coef> 3.37562e-12 </coef> *T^4+
              <coef> 992472 </coef> /T)
            </IVT_Tp> 
            <limit_Tp> T = {1500, 6000} p = {0, inf} </limit_Tp>
          </compound_Tp> 
        </species>
        <func_x> +x(Bi4)* </func_x>
        <species> Bi4
          <compound_Tp class=func_Tp>
            <IVT_Tp class=func_Tp>
              (
              <coef> 212114 </coef> +
              <coef> 160.706 </coef> *T+
              <coef> -82.8812 </coef> *T*log(T)+
              <coef> -0.000409382 </coef> *T^2+
              <coef> 1.50816e-07 </coef> *T^3+
              <coef> -2.80655e-11 </coef> *T^4+
              <coef> 53658.2 </coef> /T)
            </IVT_Tp> 
            <limit_Tp> T = {298.15, 1500} p = {0, inf} </limit_Tp>
            <IVT_Tp class=func_Tp>
              (
              <coef> 279023 </coef> +
              <coef> -109.731 </coef> *T+
              <coef> -50.0227 </coef> *T*log(T)+
              <coef> -0.00285565 </coef> *T^2+
              <coef> -6.87913e-07 </coef> *T^3+
              <coef> 5.04777e-11 </coef> *T^4+
              <coef> -2.05995e+07 </coef> /T)
            </IVT_Tp> 
            <limit_Tp> T = {1500, 6000} p = {0, inf} </limit_Tp>
          </compound_Tp> 
        </species>
        <func_x> +x(Se2)* </func_x>
        <species> Se2
          <compound_Tp class=func_Tp>
            <IVT_Tp class=func_Tp>
              (
              <coef> 128689 </coef> +
              <coef> 99.9651 </coef> *T+
              <coef> -52.967 </coef> *T*log(T)+
              <coef> 0.0173017 </coef> *T^2+
              <coef> -4.73855e-06 </coef> *T^3+
              <coef> 6.15124e-10 </coef> *T^4+
              <coef> 143169 </coef> /T)
            </IVT_Tp> 
            <limit_Tp> T = {298.15, 1500} p = {0, inf} </limit_Tp>
            <IVT_Tp class=func_Tp>
              (
              <coef> 130661 </coef> +
              <coef> 13.0007 </coef> *T+
              <coef> -38.7758 </coef> *T*log(T)+
              <coef> -0.000819316 </coef> *T^2+
              <coef> 1.15058e-08 </coef> *T^3+
              <coef> 1.86849e-12 </coef> *T^4+
              <coef> 1.10946e+06 </coef> /T)
            </IVT_Tp> 
            <limit_Tp> T = {1500, 6000} p = {0, inf} </limit_Tp>
          </compound_Tp> 
        </species>
        <func_x> +x(Se3)* </func_x>
        <species> Se3
          <compound_Tp class=func_Tp>
            <IVT_Tp class=func_Tp>
              (
              <coef> 159789 </coef> +
              <coef> 74.3728 </coef> *T+
              <coef> -57.288 </coef> *T*log(T)+
              <coef> -0.00098465 </coef> *T^2+
              <coef> 2.57277e-07 </coef> *T^3+
              <coef> -3.46006e-11 </coef> *T^4+
              <coef> 144376 </coef> /T)
            </IVT_Tp> 
            <limit_Tp> T = {298.15, 1500} p = {0, inf} </limit_Tp>
            <IVT_Tp class=func_Tp>
              (
              <coef> 159375 </coef> +
              <coef> 80.2829 </coef> *T+
              <coef> -58.1994 </coef> *T*log(T)+
              <coef> -3.12111e-07 </coef> *T^2+
              <coef> 1.52335e-11 </coef> *T^3+
              <coef> -3.60921e-16 </coef> *T^4+
              <coef> 182383 </coef> /T)
            </IVT_Tp> 
            <limit_Tp> T = {1500, 6000} p = {0, inf} </limit_Tp>
          </compound_Tp> 
        </species>
        <func_x> +x(Se4)* </func_x>
        <species> Se4
          <compound_Tp class=func_Tp>
            <IVT_Tp class=func_Tp>
              (
              <coef> 139606 </coef> +
              <coef> 205.397 </coef> *T+
              <coef> -81.5429 </coef> *T*log(T)+
              <coef> -0.00172525 </coef> *T^2+
              <coef> 4.50447e-07 </coef> *T^3+
              <coef> -6.05457e-11 </coef> *T^4+
              <coef> 245066 </coef> /T)
            </IVT_Tp> 
            <limit_Tp> T = {298.15, 1500} p = {0, inf} </limit_Tp>
            <IVT_Tp class=func_Tp>
              (
              <coef> 138884 </coef> +
              <coef> 215.741 </coef> *T+
              <coef> -83.1383 </coef> *T*log(T)+
              <coef> -1.52385e-06 </coef> *T^2+
              <coef> 9.60734e-11 </coef> *T^3+
              <coef> -3.12685e-15 </coef> *T^4+
              <coef> 310607 </coef> /T)
            </IVT_Tp> 
            <limit_Tp> T = {1500, 6000} p = {0, inf} </limit_Tp>
          </compound_Tp> 
        </species>
        <func_x> +x(Se5)* </func_x>
        <species> Se5
          <compound_Tp class=func_Tp>
            <IVT_Tp class=func_Tp>
              (
              <coef> 123067 </coef> +
              <coef> 304.395 </coef> *T+
              <coef> -107.162 </coef> *T*log(T)+
              <coef> -0.00100144 </coef> *T^2+
              <coef> 2.62304e-07 </coef> *T^3+
              <coef> -3.53422e-11 </coef> *T^4+
              <coef> 228923 </coef> /T)
            </IVT_Tp> 
            <limit_Tp> T = {298.15, 1500} p = {0, inf} </limit_Tp>
            <IVT_Tp class=func_Tp>
              (
              <coef> 122653 </coef> +
              <coef> 310.36 </coef> *T+
              <coef> -108.083 </coef> *T*log(T)+
              <coef> -1.24385e-06 </coef> *T^2+
              <coef> 8.69881e-11 </coef> *T^3+
              <coef> -3.2108e-15 </coef> *T^4+
              <coef> 265847 </coef> /T)
            </IVT_Tp> 
            <limit_Tp> T = {1500, 6000} p = {0, inf} </limit_Tp>
          </compound_Tp> 
        </species>
        <func_x> +x(Se6)* </func_x>
        <species> Se6
          <compound_Tp class=func_Tp>
            <IVT_Tp class=func_Tp>
              (
              <coef> 96728.7 </coef> +
              <coef> 452.489 </coef> *T+
              <coef> -132.249 </coef> *T*log(T)+
              <coef> -0.000846196 </coef> *T^2+
              <coef> 2.21597e-07 </coef> *T^3+
              <coef> -2.98465e-11 </coef> *T^4+
              <coef> 265562 </coef> /T)
            </IVT_Tp> 
            <limit_Tp> T = {298.15, 1500} p = {0, inf} </limit_Tp>
            <IVT_Tp class=func_Tp>
              (
              <coef> 96386.9 </coef> +
              <coef> 457.488 </coef> *T+
              <coef> -133.022 </coef> *T*log(T)+
              <coef> -2.87329e-06 </coef> *T^2+
              <coef> 2.14953e-10 </coef> *T^3+
              <coef> -8.17164e-15 </coef> *T^4+
              <coef> 294936 </coef> /T)
            </IVT_Tp> 
            <limit_Tp> T = {1500, 6000} p = {0, inf} </limit_Tp>
          </compound_Tp> 
        </species>
        <func_x> +x(Se7)* </func_x>
        <species> Se7
          <compound_Tp class=func_Tp>
            <IVT_Tp class=func_Tp>
              (
              <coef> 99421.2 </coef> +
              <coef> 561.666 </coef> *T+
              <coef> -156.927 </coef> *T*log(T)+
              <coef> -0.00113367 </coef> *T^2+
              <coef> 2.96981e-07 </coef> *T^3+
              <coef> -4.00047e-11 </coef> *T^4+
              <coef> 328173 </coef> /T)
            </IVT_Tp> 
            <limit_Tp> T = {298.15, 1500} p = {0, inf} </limit_Tp>
            <IVT_Tp class=func_Tp>
              (
              <coef> 98952.4 </coef> +
              <coef> 568.42 </coef> *T+
              <coef> -157.97 </coef> *T*log(T)+
              <coef> -1.2193e-06 </coef> *T^2+
              <coef> 8.2395e-11 </coef> *T^3+
              <coef> -2.87662e-15 </coef> *T^4+
              <coef> 370452 </coef> /T)
            </IVT_Tp> 
            <limit_Tp> T = {1500, 6000} p = {0, inf} </limit_Tp>
          </compound_Tp> 
        </species>
        <func_x> +x(Se8)* </func_x>
        <species> Se8
          <compound_Tp class=func_Tp>
            <IVT_Tp class=func_Tp>
              (
              <coef> 100885 </coef> +
              <coef> 688.655 </coef> *T+
              <coef> -181.634 </coef> *T*log(T)+
              <coef> -0.00138829 </coef> *T^2+
              <coef> 3.6337e-07 </coef> *T^3+
              <coef> -4.89099e-11 </coef> *T^4+
              <coef> 378455 </coef> /T)
            </IVT_Tp> 
            <limit_Tp> T = {298.15, 1500} p = {0, inf} </limit_Tp>
            <IVT_Tp class=func_Tp>
              (
              <coef> 100296 </coef> +
              <coef> 697.008 </coef> *T+
              <coef> -182.922 </coef> *T*log(T)+
              <coef> 1.74412e-06 </coef> *T^2+
              <coef> -1.53018e-10 </coef> *T^3+
              <coef> 6.39937e-15 </coef> *T^4+
              <coef> 433755 </coef> /T)
            </IVT_Tp> 
            <limit_Tp> T = {1500, 6000} p = {0, inf} </limit_Tp>
          </compound_Tp> 
        </species>
        <func_x> +x(BiSe)* </func_x>
        <species> BiSe
          <compound_Tp class=func_Tp>
            <IVT_Tp class=func_Tp>
              (
              <coef> 162470 </coef> +
              <coef> -20.4452 </coef> *T+
              <coef> -37.0299 </coef> *T*log(T)+
              <coef> -0.000755913 </coef> *T^2+
              <coef> 1.86633e-07 </coef> *T^3+
              <coef> -3.72122e-11 </coef> *T^4+
              <coef> 40545.7 </coef> /T)
            </IVT_Tp> 
            <limit_Tp> T = {298.15, 1500} p = {0, inf} </limit_Tp>
            <IVT_Tp class=func_Tp>
              (
              <coef> 169206 </coef> +
              <coef> -29.1018 </coef> *T+
              <coef> -36.7002 </coef> *T*log(T)+
              <coef> 0.00207156 </coef> *T^2+
              <coef> -4.94165e-07 </coef> *T^3+
              <coef> 2.60853e-11 </coef> *T^4+
              <coef> -2.5877e+06 </coef> /T)
            </IVT_Tp> 
            <limit_Tp> T = {1500, 6000} p = {0, inf} </limit_Tp>
          </compound_Tp> 
        </species>
      </Reference> 
      <IdealMixing class=FuncTpx>
        +R*T*x(Bi)*log(x(Bi))
        +R*T*x(Se)*log(x(Se))
        +R*T*x(Bi2)*log(x(Bi2))
        +R*T*x(Bi3)*log(x(Bi3))
        +R*T*x(Bi4)*log(x(Bi4))
        +R*T*x(Se2)*log(x(Se2))
        +R*T*x(Se3)*log(x(Se3))
        +R*T*x(Se4)*log(x(Se4))
        +R*T*x(Se5)*log(x(Se5))
        +R*T*x(Se6)*log(x(Se6))
        +R*T*x(Se7)*log(x(Se7))
        +R*T*x(Se8)*log(x(Se8))
        +R*T*x(BiSe)*log(x(BiSe))
      </IdealMixing> 
    </SimpleSolution> 
</AssociatedSolution> 

