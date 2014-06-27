<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<PointPhase class=phase id=BaCu13_s>
  <species> Ba0.07143Cu0.92857
    <Cp_zero class=func_Tp>
      (
      <coef id=s113a unknown=1> -924.6 </coef> +
      <coef id=s113b unknown=0> 0 </coef> *T)
    </Cp_zero> 
    <ref_plane>
      0.07143
      <species id=Ba_s> Ba
        <compound_Tp class=func_Tp>
          <SGTE_Tp class=func_Tp>
            (
            <coef> -17685.226 </coef> +
            <coef> 233.78606 </coef> *T+
            <coef> -42.889 </coef> *T*log(T)+
            <coef> -0.0018314 </coef> *T^2+
            <coef> -9.5e-11 </coef> *T^3+
            <coef> 0 </coef> *T^7+
            <coef> 705880 </coef> /T+
            <coef> 0 </coef> /T^9)
          </SGTE_Tp> 
          <limit_Tp> T = {298.15, 1000} p = {0, inf} </limit_Tp>
          <SGTE_Tp class=func_Tp>
            (
            <coef> -64873.614 </coef> +
            <coef> 608.188389 </coef> *T+
            <coef> -94.2824199 </coef> *T*log(T)+
            <coef> 0.019504772 </coef> *T^2+
            <coef> -1.051353e-06 </coef> *T^3+
            <coef> 0 </coef> *T^7+
            <coef> 8220192 </coef> /T+
            <coef> 0 </coef> /T^9)
          </SGTE_Tp> 
          <limit_Tp> T = {1000, 2995} p = {0, inf} </limit_Tp>
        </compound_Tp> 
      </species>
      0.92857
      <species id=Cu_s> Cu
        <compound_Tp class=func_Tp>
          <SGTE_Tp class=func_Tp>
            (
            <coef> -7770.458 </coef> +
            <coef> 130.485235 </coef> *T+
            <coef> -24.112392 </coef> *T*log(T)+
            <coef> -0.00265684 </coef> *T^2+
            <coef> 1.29223e-07 </coef> *T^3+
            <coef> 0 </coef> *T^7+
            <coef> 52478 </coef> /T+
            <coef> 0 </coef> /T^9)
          </SGTE_Tp> 
          <limit_Tp> T = {298.15, 1357.77} p = {0, inf} </limit_Tp>
          <SGTE_Tp class=func_Tp>
            (
            <coef> -13542.026 </coef> +
            <coef> 183.803828 </coef> *T+
            <coef> -31.38 </coef> *T*log(T)+
            <coef> 0 </coef> *T^2+
            <coef> 0 </coef> *T^3+
            <coef> 0 </coef> *T^7+
            <coef> 0 </coef> /T+
            <coef> 3.642e+29 </coef> /T^9)
          </SGTE_Tp> 
          <limit_Tp> T = {1357.77, 3200} p = {0, inf} </limit_Tp>
        </compound_Tp> 
      </species>
    </ref_plane>
  </species>
</PointPhase> 
<PointPhase class=phase id=BaCu_s>
  <species> Ba0.5Cu0.5
    <Cp_zero class=func_Tp>
      (
      <coef id=s11a unknown=1> -2594 </coef> +
      <coef id=s11b unknown=0> 0 </coef> *T)
    </Cp_zero> 
    <ref_plane>
      0.5
      <species IDREF=Ba_s></species>
      0.5
      <species IDREF=Cu_s></species>
    </ref_plane>
  </species>
</PointPhase> 
<PointPhase class=phase id=Ba_s>
  <species IDREF=Ba_s></species>
</PointPhase> 
<PointPhase class=phase id=Cu_s>
  <species IDREF=Cu_s></species>
</PointPhase> 
<SimpleSolution class=phase id=L>
  <components> Ba Cu </components>
  <Reference class=FuncTpx>
    <func_x> +x(Ba)* </func_x>
    <species id=Ba_l> Ba
      <compound_Tp class=func_Tp>
        <SGTE_Tp class=func_Tp>
          (
          <coef> -9738.988 </coef> +
          <coef> 229.540143 </coef> *T+
          <coef> -43.4961089 </coef> *T*log(T)+
          <coef> -0.002346416 </coef> *T^2+
          <coef> 9.91223e-07 </coef> *T^3+
          <coef> 0 </coef> *T^7+
          <coef> 723016 </coef> /T+
          <coef> 0 </coef> /T^9)
        </SGTE_Tp> 
        <limit_Tp> T = {298.15, 1000} p = {0, inf} </limit_Tp>
        <SGTE_Tp class=func_Tp>
          (
          <coef> -7381.093 </coef> +
          <coef> 235.49642 </coef> *T+
          <coef> -45.103 </coef> *T*log(T)+
          <coef> 0.002154 </coef> *T^2+
          <coef> 2.7e-11 </coef> *T^3+
          <coef> 0 </coef> *T^7+
          <coef> -365 </coef> /T+
          <coef> 0 </coef> /T^9)
        </SGTE_Tp> 
        <limit_Tp> T = {1000, 2995} p = {0, inf} </limit_Tp>
      </compound_Tp> 
    </species>
    <func_x> +x(Cu)* </func_x>
    <species id=Cu_l> Cu
      <compound_Tp class=func_Tp>
        <SGTE_Tp class=func_Tp>
          (
          <coef> 5194.277 </coef> +
          <coef> 120.973331 </coef> *T+
          <coef> -24.112392 </coef> *T*log(T)+
          <coef> -0.00265684 </coef> *T^2+
          <coef> 1.29223e-07 </coef> *T^3+
          <coef> -5.849e-21 </coef> *T^7+
          <coef> 52478 </coef> /T+
          <coef> 0 </coef> /T^9)
        </SGTE_Tp> 
        <limit_Tp> T = {298.15, 1357.77} p = {0, inf} </limit_Tp>
        <SGTE_Tp class=func_Tp>
          (
          <coef> -46.545 </coef> +
          <coef> 173.881484 </coef> *T+
          <coef> -31.38 </coef> *T*log(T)+
          <coef> 0 </coef> *T^2+
          <coef> 0 </coef> *T^3+
          <coef> 0 </coef> *T^7+
          <coef> 0 </coef> /T+
          <coef> 0 </coef> /T^9)
        </SGTE_Tp> 
        <limit_Tp> T = {1357.77, 3200} p = {0, inf} </limit_Tp>
      </compound_Tp> 
    </species>
  </Reference> 
  <IdealMixing class=FuncTpx>
    +R*T*x(Ba)*log(x(Ba))
    +R*T*x(Cu)*log(x(Cu))
  </IdealMixing> 
  <RedlichKister class=FuncTpx formalism=NoChange>
    +x(Ba)*x(Cu)*(
      <func_x> + </func_x> 
      <Cp_zero class=func_Tp>
        (
        <coef id=L0a unknown=1> -7191 </coef> +
        <coef id=L0b unknown=1> 3.363 </coef> *T)
      </Cp_zero> 
      <func_x> +(x(Ba)-x(Cu))* </func_x> 
      <Cp_zero class=func_Tp>
        (
        <coef id=L1a unknown=0> 0 </coef> +
        <coef id=L1b unknown=0> 0 </coef> *T)
      </Cp_zero> 
      <func_x> +(x(Ba)-x(Cu))^2* </func_x> 
      <Cp_zero class=func_Tp>
        (
        <coef id=L2a unknown=0> 0 </coef> +
        <coef id=L2b unknown=0> 0 </coef> *T)
      </Cp_zero> 
    )
  </RedlichKister> 
</SimpleSolution> 



