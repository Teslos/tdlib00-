<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<associated_solution class=phase id=L
  debug=0>
  <components> Bi Se </components>
  <Reference class=FuncTpx>
    <func_x> +x(Bi)* </func_x>
    <species> Bi
      <compound_Tp class=func_Tp>
        <SGTE_Tp class=func_Tp>
          (
          <coef> 3428.29 </coef> +
          <coef> 107.782 </coef> *T+
          <coef> -28.4097 </coef> *T*log(T)+
          <coef> 0.0123389 </coef> *T^2+
          <coef> -8.3816e-06 </coef> *T^3+
          <coef> -5.995e-19 </coef> *T^7+
          <coef> 0 </coef> /T+
          <coef> 0 </coef> /T^9)
        </SGTE_Tp> 
        <limit_Tp> T = {250, 544.55} p = {0, inf} </limit_Tp>
        <SGTE_Tp class=func_Tp>
          (
          <coef> 41544.3 </coef> +
          <coef> -414.461 </coef> *T+
          <coef> 51.8557 </coef> *T*log(T)+
          <coef> -0.0753112 </coef> *T^2+
          <coef> 1.34999e-05 </coef> *T^3+
          <coef> 0 </coef> *T^7+
          <coef> -3.61617e+06 </coef> /T+
          <coef> 0 </coef> /T^9)
        </SGTE_Tp> 
        <limit_Tp> T = {544.55, 800} p = {0, inf} </limit_Tp>
        <SGTE_Tp class=func_Tp>
          (
          <coef> 290.595 </coef> +
          <coef> 161.739 </coef> *T+
          <coef> -35.9824 </coef> *T*log(T)+
          <coef> 0.0074266 </coef> *T^2+
          <coef> -1.046e-06 </coef> *T^3+
          <coef> 0 </coef> *T^7+
          <coef> 0 </coef> /T+
          <coef> 0 </coef> /T^9)
        </SGTE_Tp> 
        <limit_Tp> T = {800, 1200} p = {0, inf} </limit_Tp>
        <SGTE_Tp class=func_Tp>
          (
          <coef> 3754.95 </coef> +
          <coef> 103.961 </coef> *T+
          <coef> -27.196 </coef> *T*log(T)+
          <coef> 0 </coef> *T^2+
          <coef> 0 </coef> *T^3+
          <coef> 0 </coef> *T^7+
          <coef> 0 </coef> /T+
          <coef> 0 </coef> /T^9)
        </SGTE_Tp> 
        <limit_Tp> T = {1200, 3000} p = {0, inf} </limit_Tp>
      </compound_Tp> 
    </species>
    <func_x> +x(Se)* </func_x>
    <species> Se
      <compound_Tp class=func_Tp>
        <SGTE_Tp class=func_Tp>
          (
          <coef> 50533.3 </coef> +
          <coef> -1178.29 </coef> *T+
          <coef> 194.107 </coef> *T*log(T)+
          <coef> -0.390269 </coef> *T^2+
          <coef> 0.000119219 </coef> *T^3+
          <coef> 0 </coef> *T^7+
          <coef> -2.2244e+06 </coef> /T+
          <coef> 0 </coef> /T^9)
        </SGTE_Tp> 
        <limit_Tp> T = {250, 494} p = {0, inf} </limit_Tp>
        <SGTE_Tp class=func_Tp>
          (
          <coef> -5228.3 </coef> +
          <coef> 183.726 </coef> *T+
          <coef> -35.1456 </coef> *T*log(T)+
          <coef> 0 </coef> *T^2+
          <coef> 0 </coef> *T^3+
          <coef> 0 </coef> *T^7+
          <coef> 0 </coef> /T+
          <coef> 0 </coef> /T^9)
        </SGTE_Tp> 
        <limit_Tp> T = {494, 1300} p = {0, inf} </limit_Tp>
      </compound_Tp> 
    </species>
  </Reference> 
  <internal_solution>
  <SimpleSolution class=phase>
    <components> Bi Se Bi0.4Se0.6 </components>
    <Reference class=FuncTpx>
      <func_x> +x(Bi)* </func_x>
      <species> Bi
      </species>
      <func_x> +x(Se)* </func_x>
      <species> Se
      </species>
      <func_x> +x(Bi0.4Se0.6)* </func_x>
      <species> Bi0.4Se0.6
        <Cp_const class=func_Tp>
          (
          <coef> -23849.1 </coef> +
          <coef> -0.951154 </coef> *T+
          <coef> 1.01578 </coef> *T*log(T))
        </Cp_const> 
      </species>
    </Reference> 
    <IdealMixing class=FuncTpx>
      +R*T*x(Bi)*log(x(Bi))
      +R*T*x(Se)*log(x(Se))
      +R*T*x(Bi0.4Se0.6)*log(x(Bi0.4Se0.6))
    </IdealMixing> 
    <RedlichKister class=FuncTpx formalism=Muggianu>
      +x(Bi)*x(Bi0.4Se0.6)*(
        <func_x> + </func_x> 
        <Cp_const class=func_Tp>
          (
          <coef> 1154.49 </coef> +
          <coef> 2.47494 </coef> *T+
          <coef> 0 </coef> *T*log(T))
        </Cp_const> 
        <func_x> +(x(Bi)-x(Bi0.4Se0.6))* </func_x> 
        <Cp_const class=func_Tp>
          (
          <coef> 11298.6 </coef> +
          <coef> -11.9604 </coef> *T+
          <coef> 0 </coef> *T*log(T))
        </Cp_const> 
      )
    </RedlichKister> 
    <RedlichKister class=FuncTpx formalism=Muggianu>
      +x(Se)*x(Bi0.4Se0.6)*(
        <func_x> + </func_x> 
        <Cp_const class=func_Tp>
          (
          <coef> 15590.5 </coef> +
          <coef> 2.34574 </coef> *T+
          <coef> 0 </coef> *T*log(T))
        </Cp_const> 
        <func_x> +(x(Se)-x(Bi0.4Se0.6))* </func_x> 
        <null_Tp class=func_Tp>
          (0)
        </null_Tp> 
        <func_x> +(x(Se)-x(Bi0.4Se0.6))^2* </func_x> 
        <Cp_const class=func_Tp>
          (
          <coef> -21418.4 </coef> +
          <coef> 21.1963 </coef> *T+
          <coef> 0 </coef> *T*log(T))
        </Cp_const> 
      )
    </RedlichKister> 
  </SimpleSolution> 
  </internal_solution>
</associated_solution> 
<PointPhase class=phase id=s01>
  <species id=Se_s> Se
    <compound_Tp class=func_Tp>
      <SGTE_Tp class=func_Tp>
        (
        <coef> -9376.37 </coef> +
        <coef> 174.206 </coef> *T+
        <coef> -33.6527 </coef> *T*log(T)+
        <coef> 0.0242431 </coef> *T^2+
        <coef> -1.53185e-05 </coef> *T^3+
        <coef> 0 </coef> *T^7+
        <coef> 102249 </coef> /T+
        <coef> 0 </coef> /T^9)
      </SGTE_Tp> 
      <limit_Tp> T = {250, 494} p = {0, inf} </limit_Tp>
      <SGTE_Tp class=func_Tp>
        (
        <coef> -37546.1 </coef> +
        <coef> 507.112 </coef> *T+
        <coef> -81.2007 </coef> *T*log(T)+
        <coef> 0.0371449 </coef> *T^2+
        <coef> -5.61103e-06 </coef> *T^3+
        <coef> 0 </coef> *T^7+
        <coef> 2.61426e+06 </coef> /T+
        <coef> 0 </coef> /T^9)
      </SGTE_Tp> 
      <limit_Tp> T = {494, 800} p = {0, inf} </limit_Tp>
      <SGTE_Tp class=func_Tp>
        (
        <coef> -12193.1 </coef> +
        <coef> 197.77 </coef> *T+
        <coef> -35.1456 </coef> *T*log(T)+
        <coef> 0 </coef> *T^2+
        <coef> 0 </coef> *T^3+
        <coef> 0 </coef> *T^7+
        <coef> 0 </coef> /T+
        <coef> 0 </coef> /T^9)
      </SGTE_Tp> 
      <limit_Tp> T = {800, 1300} p = {0, inf} </limit_Tp>
    </compound_Tp> 
  </species>
</PointPhase> 
<PointPhase class=phase id=s10>
  <species id=Bi_s> Bi
    <compound_Tp class=func_Tp>
      <SGTE_Tp class=func_Tp>
        (
        <coef> -7817.78 </coef> +
        <coef> 128.419 </coef> *T+
        <coef> -28.4097 </coef> *T*log(T)+
        <coef> 0.0123389 </coef> *T^2+
        <coef> -8.3816e-06 </coef> *T^3+
        <coef> 0 </coef> *T^7+
        <coef> 0 </coef> /T+
        <coef> 0 </coef> /T^9)
      </SGTE_Tp> 
      <limit_Tp> T = {250, 544.55} p = {0, inf} </limit_Tp>
      <SGTE_Tp class=func_Tp>
        (
        <coef> 30208 </coef> +
        <coef> -393.65 </coef> *T+
        <coef> 51.8557 </coef> *T*log(T)+
        <coef> -0.0753112 </coef> *T^2+
        <coef> 1.34999e-05 </coef> *T^3+
        <coef> 0 </coef> *T^7+
        <coef> -3.61617e+06 </coef> /T+
        <coef> 1.661e+25 </coef> /T^9)
      </SGTE_Tp> 
      <limit_Tp> T = {544.55, 800} p = {0, inf} </limit_Tp>
      <SGTE_Tp class=func_Tp>
        (
        <coef> -11045.7 </coef> +
        <coef> 182.549 </coef> *T+
        <coef> -35.9824 </coef> *T*log(T)+
        <coef> 0.0074266 </coef> *T^2+
        <coef> -1.046e-06 </coef> *T^3+
        <coef> 0 </coef> *T^7+
        <coef> 0 </coef> /T+
        <coef> 1.661e+25 </coef> /T^9)
      </SGTE_Tp> 
      <limit_Tp> T = {800, 1200} p = {0, inf} </limit_Tp>
      <SGTE_Tp class=func_Tp>
        (
        <coef> -7581.31 </coef> +
        <coef> 124.771 </coef> *T+
        <coef> -27.196 </coef> *T*log(T)+
        <coef> 0 </coef> *T^2+
        <coef> 0 </coef> *T^3+
        <coef> 0 </coef> *T^7+
        <coef> 0 </coef> /T+
        <coef> 1.661e+25 </coef> /T^9)
      </SGTE_Tp> 
      <limit_Tp> T = {1200, 3000} p = {0, inf} </limit_Tp>
    </compound_Tp> 
  </species>
</PointPhase> 
<SimpleSolution class=phase id=s11>
  <components> Bi Se </components>
  <Reference class=FuncTpx>
    <func_x> +x(Bi)* </func_x>
    <species> Bi
      <Cp_BB2_Tref class=func_Tp
        To=300>
        (
        <coef> 27513.8 </coef> -
        <coef> 14.8122 </coef> *T+
        <coef id=s11_Cpa> 4.01738 </coef> *(T-To-T*log(T/To))+4*
        <coef id=s11_Cpb> -34.1601 </coef> *(sqrt(T)-sqrt(To)))
      </Cp_BB2_Tref> 
      <ref_plane>
        1
        <species id=Bi_ref> Bi_ref
          <Cp_BB2_Tref class=func_Tp
            To=300>
            (
            <coef> 47.2354 </coef> -
            <coef> 56.8931 </coef> *T+
            <coef> 25.5324 </coef> *(T-To-T*log(T/To))+4*
            <coef> 0 </coef> *(sqrt(T)-sqrt(To)))
          </Cp_BB2_Tref> 
        </species>
      </ref_plane>
    </species>
    <func_x> +x(Se)* </func_x>
    <species> Se
      <Cp_BB2_Tref class=func_Tp
        To=300>
        (
        <coef> -17167.2 </coef> -
        <coef> -13.3682 </coef> *T+
        <coef IDREF=s11_Cpa></coef> *(T-To-T*log(T/To))+4*
        <coef IDREF=s11_Cpb></coef> *(sqrt(T)-sqrt(To)))
      </Cp_BB2_Tref> 
      <ref_plane>
        1
        <species id=Se_ref> Se_ref
          <Cp_BB2_Tref class=func_Tp
            To=300>
            (
            <coef> 46.42 </coef> -
            <coef> 42.1206 </coef> *T+
            <coef> 25.1066 </coef> *(T-To-T*log(T/To))+4*
            <coef> 0 </coef> *(sqrt(T)-sqrt(To)))
          </Cp_BB2_Tref> 
        </species>
      </ref_plane>
    </species>
  </Reference> 
  <IdealMixing class=FuncTpx>
    +R*T*x(Bi)*log(x(Bi))
    +R*T*x(Se)*log(x(Se))
  </IdealMixing> 
  <RedlichKister class=FuncTpx formalism=NoChange>
    +x(Bi)*x(Se)*(
      <func_x> + </func_x> 
      <Cp_BB2 class=func_Tp>
        (
        <coef> -124147 </coef> +
        <coef> 41.6101 </coef> *T+
        <coef> 0 </coef> *T*log(T)+
        <coef> 0 </coef> *sqrt(T))
      </Cp_BB2> 
    )
  </RedlichKister> 
</SimpleSolution> 
<PointPhase class=phase id=s23>
  <species> Bi0.4Se0.6
    <Cp_BB2_Tref class=func_Tp
      To=300>
      (
      <coef> -31709.4 </coef> -
      <coef> -14.3682 </coef> *T+
      <coef> 7.30398 </coef> *(T-To-T*log(T/To))+4*
      <coef> -107.547 </coef> *(sqrt(T)-sqrt(To)))
    </Cp_BB2_Tref> 
    <ref_plane>
      0.4
      <species IDREF=Bi_ref></species>
      0.6
      <species IDREF=Se_ref></species>
    </ref_plane>
  </species>
</PointPhase> 
<PointPhase class=phase id=s32>
  <species> Bi0.6Se0.4
    <Cp_BB2_Tref class=func_Tp
      To=300>
      (
      <coef> -18422.6 </coef> -
      <coef> 4.65502 </coef> *T+
      <coef> 0 </coef> *(T-To-T*log(T/To))+4*
      <coef> 0 </coef> *(sqrt(T)-sqrt(To)))
    </Cp_BB2_Tref> 
    <ref_plane>
      0.6
      <species IDREF=Bi_ref></species>
      0.4
      <species IDREF=Se_ref></species>
    </ref_plane>
  </species>
</PointPhase> 


