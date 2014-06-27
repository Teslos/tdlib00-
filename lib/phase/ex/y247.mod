<CuOx_ordered_plane class=phase id=Y247>
  <components> Y2Ba4Cu7O14 Y2Ba4Cu7O15 </components>
  <Reference class=FuncTpx>
    <func_x> +x(Y2Ba4Cu7O14)* </func_x>
    <species id=Y2Ba4Cu7O14> Y2Ba4Cu7O14
      <ref_plane>
        1
        <species id=Y2O3> Y2O3
          <Cp_BB4 class=func_Tp>
            (
            <coef id=a unknown=0> -30047.8 </coef> +
            <coef id=b unknown=0> 954.105 </coef> *T+
            <coef id=c unknown=0> -146.996 </coef> *T*log(T)+
            <coef id=d unknown=0> -2120.1 </coef> *sqrt(T)+
            <coef id=e unknown=0> 737146 </coef> /T+
            <coef id=f unknown=0> -1.24419e+07 </coef> /T^2)
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
            <coef> -9.96343e+06 </coef> /T^2)
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
    <species id=Y2Ba4Cu7O15> Y2Ba4Cu7O15
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
            <coef> -1.3651e+07 </coef> /T^2)
          </Cp_BB4> 
        </species>
      </ref_plane>
    </species>
  </Reference> 
  2*R*T*(z*log(z)+(1-z)*log(1-z))
  <func_x> + </func_x> 
    <Cp_const class=func_Tp>
      (
      <coef id=2A1 unknown=1> -147556 </coef> +
      <coef id=2B1 unknown=1> -10.4806 </coef> *T+
      <coef id=2C1 unknown=0> 0 </coef> *T*log(T))
    </Cp_const> 
  <func_x> +z* </func_x> 
  <Cp_BB4 class=func_Tp>
    (
    <coef id=2A2 unknown=1> -74531 </coef> +
    <coef id=2B2 unknown=1> 1056.5 </coef> *T+
    <coef id=2C2 unknown=1> -107.411 </coef> *T*log(T)+
    <coef id=2D2 unknown=1> -7908.52 </coef> *sqrt(T)+
    <coef> 0 </coef> /T+
    <coef> 0 </coef> /T^2)
  </Cp_BB4> 
  <func_x> +z*(1-z)^1* </func_x> 
  <Cp_const class=func_Tp>
    (
    <coef id=2ah1 unknown=1> -3725.4 </coef> +
    <coef id=2as1 unknown=1> 12.4786 </coef> *T+
    <coef id=2ac1 unknown=0> 0 </coef> *T*log(T))
  </Cp_const> 
  <func_x> +z*(1-z)^2* </func_x> 
  <Cp_const class=func_Tp>
    (
    <coef id=2ah2 unknown=1> -43801.5 </coef> +
    <coef id=2as2 unknown=1> 39.792 </coef> *T+
    <coef id=2ac2 unknown=0> 0 </coef> *T*log(T))
  </Cp_const> 
  <func_x> +z*(1-z)^3* </func_x> 
  <Cp_const class=func_Tp>
    (
    <coef id=2ah3 unknown=0> 0 </coef> +
    <coef id=2as3 unknown=0> 0 </coef> *T+
    <coef id=2ac3 unknown=0> 0 </coef> *T*log(T))
  </Cp_const> 
  <func_x> +z*(1-z)^4* </func_x> 
  <Cp_const class=func_Tp>
    (
    <coef id=2ah4 unknown=0> 0 </coef> +
    <coef id=2as4 unknown=0> 0 </coef> *T+
    <coef id=2ac4 unknown=0> 0 </coef> *T*log(T))
  </Cp_const> 
</CuOx_ordered_plane> 

