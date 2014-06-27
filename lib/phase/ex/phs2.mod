<pol_solution class=phase id=L>
  <components> 
    <species> Ba_l
      <compound_Tp class=func_Tp> 
        <SGTE_Tp class=func_Tp> ( 
          <v name=a unknown=0>-9738.988    </v> +
          <v name=b unknown=0>+229.540143  </v> *T+
          <v name=c unknown=0>-43.4961089  </v> *T*log(T)+
          <v name=d unknown=0>-2.346416e-3 </v> *T^2+
          <v name=e unknown=0>+0.991223e-6 </v> *T^3+
          <v name=f unknown=0>0            </v> *T^7+
          <v name=g unknown=0>+723016      </v> /T+
          <v name=h unknown=0> 0 </v> /T^9)
        </SGTE_Tp> 
        <limit_Tp> T = {298.15, 1000} p = {0, inf} </limit_Tp>
        <SGTE_Tp class=func_Tp> ( 
          <v name=a unknown=0> -7381.093    </v> +
          <v name=b unknown=0> +235.49642   </v> *T+
          <v name=c unknown=0> -45.103      </v> *T*log(T)+
          <v name=d unknown=0> +2.154e-3    </v> *T^2+
          <v name=e unknown=0> +0.000027e-6 </v> *T^3+
          <v name=f unknown=0> 0            </v> *T^7+
          <v name=g unknown=0> -365         </v> /T+
          <v name=h unknown=0> 0 </v> /T^9)
        </SGTE_Tp> 
        <limit_Tp> T = {1000, 2995} p = {0, inf} </limit_Tp>
      </compound_Tp>
    </species>
    <species> Cu_l
      <compound_Tp class=func_Tp> 
        <SGTE_Tp class=func_Tp> ( 
          <v name=a unknown=0> 5194.277    </v> +
          <v name=b unknown=0> +120.973331 </v> *T+
          <v name=c unknown=0> -24.112392  </v> *T*log(T)+
          <v name=d unknown=0> -2.65684e-3 </v> *T^2+
          <v name=e unknown=0> +0.129223e-6</v> *T^3+
          <v name=f unknown=0> -5.849e-21  </v> *T^7+
          <v name=g unknown=0> +52478      </v> /T+
          <v name=h unknown=0> 0 </v> /T^9)
        </SGTE_Tp> 
        <limit_Tp> T = {298.15, 1357.77} p = {0, inf} </limit_Tp>
        <SGTE_Tp class=func_Tp> ( 
          <v name=a unknown=0> -46.545    </v> +
          <v name=b unknown=0> +173.881484</v> *T+
          <v name=c unknown=0> -31.38     </v> *T*log(T)+
          <v name=d unknown=0> 0 </v> *T^2+
          <v name=e unknown=0> 0 </v> *T^3+
          <v name=f unknown=0> 0 </v> *T^7+
          <v name=g unknown=0> 0 </v> /T+
          <v name=h unknown=0> 0 </v> /T^9)
        </SGTE_Tp> 
        <limit_Tp> T = {1357.77, 3200} p = {0, inf} </limit_Tp>
      </compound_Tp>
    </species>
  </components>
  Gid 
  <Redlich_Kister class=pol_interaction >
    +x1*x2*(
    <func_x> +1*</func_x> 
    <Cp_zero class=func_Tp> ( 
      <v name=La0 unknown=1> -7191. </v> +
      <v name=Lb0 unknown=1> 3.363 </v> *T)
    </Cp_zero> 
    <func_x> +(x1-x2)*</func_x> 
    <Cp_zero class=func_Tp> ( 
      <v name=La1 unknown=0> 0 </v> +
      <v name=Lb1 unknown=0> 0 </v> *T)
    </Cp_zero> 
    <func_x> +(x1-x2)^2*</func_x> 
    <Cp_zero class=func_Tp> ( 
      <v name=La2 unknown=0> 0 </v> +
      <v name=Lb2 unknown=0> 0 </v> *T)
    </Cp_zero> 
  ) </Redlich_Kister>
</pol_solution>

