<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<CuOx_plane class=phase id=Y123
  debug=0>
  <components> YBa2Cu3O6 YBa2Cu3O7 </components>
  <Reference class=FuncTpx>
    <func_x> +x(YBa2Cu3O6)* </func_x>
    <species> YBa2Cu3O6
      <ref_plane>
        0.5
        <species id=Y2O3> Y2O3
          <Cp_BB4 class=func_Tp>
            (
            <coef> -30047.8 </coef> +
            <coef> 954.105 </coef> *T+
            <coef> -146.996 </coef> *T*log(T)+
            <coef> -2120.1 </coef> *sqrt(T)+
            <coef> 737146 </coef> /T+
            <coef> -1.24419e+07 </coef> /T^2)
          </Cp_BB4> 
        </species>
        2
        <species id=BaO> BaO
          <Cp_BB4 class=func_Tp>
            (
            <coef> -5092.97 </coef> +
            <coef> 463.412 </coef> *T+
            <coef> -72.028 </coef> *T*log(T)+
            <coef> -1858.56 </coef> *sqrt(T)+
            <coef> 0 </coef> /T+
            <coef> -9.96343e+06 </coef> /T^2)
          </Cp_BB4> 
        </species>
        3
        <species id=CuO> CuO
          <Cp_BB4 class=func_Tp>
            (
            <coef> -5669.13 </coef> +
            <coef> 477.648 </coef> *T+
            <coef> -69.785 </coef> *T*log(T)+
            <coef> -1801.18 </coef> *sqrt(T)+
            <coef> 61609 </coef> /T+
            <coef> 0 </coef> /T^2)
          </Cp_BB4> 
        </species>
        -0.25
        <species id=O2> O2
          <Cp_BB4 class=func_Tp>
            (
            <coef> -1776.28 </coef> +
            <coef> 132.543 </coef> *T+
            <coef> -44.978 </coef> *T*log(T)+
            <coef> -1294.17 </coef> *sqrt(T)+
            <coef> 0 </coef> /T+
            <coef> -1.3651e+07 </coef> /T^2)
          </Cp_BB4> 
        </species>
      </ref_plane>
    </species>
    <func_x> +x(YBa2Cu3O7)* </func_x>
    <species> YBa2Cu3O7
      <ref_plane>
        0.5
        <species IDREF=Y2O3></species>
        2
        <species IDREF=BaO></species>
        3
        <species IDREF=CuO></species>
        0.25
        <species IDREF=O2></species>
      </ref_plane>
    </species>
  </Reference> 
  R*T*(z*log(z)+(1-z)*log(1-z)
  +(c+x)*log(c+x)+(c-x)*log(c-x)
  +(1-c+x)*log(1-c+x)+(1-c-x)*log(1-c-x))
  <func_x> + </func_x> 
  <Cp_BB4 class=func_Tp>
    (
    <coef> -29631.8 </coef> +
    <coef> -40.8877 </coef> *T+
    <coef> 0 </coef> *T*log(T)+
    <coef> 0 </coef> *sqrt(T)+
    <coef> 0 </coef> /T+
    <coef> 0 </coef> /T^2)
  </Cp_BB4> 
  <func_x> +z* </func_x> 
  <Cp_BB4 class=func_Tp>
    (
    <coef> -86131.5 </coef> +
    <coef> 382.363 </coef> *T+
    <coef> -33.6818 </coef> *T*log(T)+
    <coef> -2097.96 </coef> *sqrt(T)+
    <coef> 0 </coef> /T+
    <coef> 0 </coef> /T^2)
  </Cp_BB4> 
  <func_x> +z*(1-z)^1* </func_x> 
  <Cp_const class=func_Tp>
    (
    <coef> 16625.8 </coef> +
    <coef> -12.8384 </coef> *T+
    <coef> 0 </coef> *T*log(T))
  </Cp_const> 
  <func_x> +z*(1-z)^2* </func_x> 
  <Cp_const class=func_Tp>
    (
    <coef> -21536.6 </coef> +
    <coef> 33.6251 </coef> *T+
    <coef> 0 </coef> *T*log(T))
  </Cp_const> 
  <func_x> +(c^2-x^2)* </func_x> 
  <Cp_zero class=func_Tp>
    (
    <coef> 5421.48 </coef> +
    <coef> 32.6017 </coef> *T)
  </Cp_zero> 
</CuOx_plane> 

