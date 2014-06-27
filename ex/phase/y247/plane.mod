<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<CuOx_ordered_plane class=phase id=Y247>
  <components> Y2Ba4Cu7O14 Y2Ba4Cu7O15 </components>
  <Reference class=FuncTpx>
    <func_x> +x(Y2Ba4Cu7O14)* </func_x>
    <species> Y2Ba4Cu7O14
      <ref_plane>
        1
        <species id=Y2O3> Y2O3
          <Cp_BB4 class=func_Tp>
            (
            <coef> -30047.8 </coef> +
            <coef> 954.105 </coef> *T+
            <coef> -146.996 </coef> *T*log(T)+
            <coef> -2120.1 </coef> *sqrt(T)+
            <coef> 737146 </coef> /T+
            <coef> -12441900 </coef> /T^2)
          </Cp_BB4> 
        </species>
        4
        <species id=BaO> BaO
          <Cp_BB4 class=func_Tp>
            (
            <coef> -5092.97 </coef> +
            <coef> 463.412 </coef> *T+
            <coef> -72.028 </coef> *T*log(T)+
            <coef> -1858.56 </coef> *sqrt(T)+
            <coef> 0 </coef> /T+
            <coef> -9963430 </coef> /T^2)
          </Cp_BB4> 
        </species>
        7
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
      </ref_plane>
    </species>
    <func_x> +x(Y2Ba4Cu7O15)* </func_x>
    <species> Y2Ba4Cu7O15
      <ref_plane>
        1
        <species IDREF=Y2O3></species>
        4
        <species IDREF=BaO></species>
        7
        <species IDREF=CuO></species>
        0.5
        <species id=O2> O2
          <Cp_BB4 class=func_Tp>
            (
            <coef> -1776.28 </coef> +
            <coef> 132.543 </coef> *T+
            <coef> -44.978 </coef> *T*log(T)+
            <coef> -1294.17 </coef> *sqrt(T)+
            <coef> 0 </coef> /T+
            <coef> -13651000 </coef> /T^2)
          </Cp_BB4> 
        </species>
      </ref_plane>
    </species>
  </Reference> 
  2*R*T*(z*log(z)+(1-z)*log(1-z))
  <func_x> + </func_x> 
    <Cp_const class=func_Tp>
      (
      <coef> -147556 </coef> +
      <coef> -10.4806 </coef> *T+
      <coef> 0 </coef> *T*log(T))
    </Cp_const> 
  <func_x> +z* </func_x> 
  <Cp_BB4 class=func_Tp>
    (
    <coef> -74531 </coef> +
    <coef> 1056.5 </coef> *T+
    <coef> -107.411 </coef> *T*log(T)+
    <coef> -7908.52 </coef> *sqrt(T)+
    <coef> 0 </coef> /T+
    <coef> 0 </coef> /T^2)
  </Cp_BB4> 
  <func_x> +z*(1-z)^1* </func_x> 
  <Cp_const class=func_Tp>
    (
    <coef> -3725.4 </coef> +
    <coef> 12.4786 </coef> *T+
    <coef> 0 </coef> *T*log(T))
  </Cp_const> 
  <func_x> +z*(1-z)^2* </func_x> 
  <Cp_const class=func_Tp>
    (
    <coef> -43801.5 </coef> +
    <coef> 39.792 </coef> *T+
    <coef> 0 </coef> *T*log(T))
  </Cp_const> 
  <func_x> +z*(1-z)^3* </func_x> 
  <Cp_const class=func_Tp>
    (
    <coef> 0 </coef> +
    <coef> 0 </coef> *T+
    <coef> 0 </coef> *T*log(T))
  </Cp_const> 
  <func_x> +z*(1-z)^4* </func_x> 
  <Cp_const class=func_Tp>
    (
    <coef> 0 </coef> +
    <coef> 0 </coef> *T+
    <coef> 0 </coef> *T*log(T))
  </Cp_const> 
</CuOx_ordered_plane> 
