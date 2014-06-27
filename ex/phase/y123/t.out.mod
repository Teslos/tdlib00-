<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<OutputFile ext=t format=file>
  <ComputeOutput class=output>
    <start>
      <var name=T value=300></var>
    </start>
    <finish>
      <var name=T value=1300 operation=GE></var>
    </finish>
    <step>
      <var name=T><convert> T+100 </convert></var>
    </step>
    <compute>
      <algorithm IDREF=twosol></algorithm>
      <input value=0> x(YBa2Cu3O7) </input>
      <input name=T> T </input>
      <output> state(all) </output>
      <output> internal(all) </output>
      <output> stable(all) </output>
      <output> G </output>
      <output> H </output>
      <output> S </output>
      <output> Cp </output>
      <output> V </output>
      <output> dVdT </output>
      <output> dVdp </output>
      <output> G(all) </output>
      <output> H(all) </output>
      <output> S(all) </output>
      <output> V(all) </output>
    </compute>
  </ComputeOutput> 
  <ComputeOutput class=output>
    <start>
      <var name=T value=300></var>
    </start>
    <finish>
      <var name=T value=1300 operation=GE></var>
    </finish>
    <step>
      <var name=T><convert> T+100 </convert></var>
    </step>
    <compute>
      <algorithm IDREF=twosol></algorithm>
      <input value=0.25> x(YBa2Cu3O7) </input>
      <input name=T> T </input>
      <output> state(all) </output>
      <output> internal(all) </output>
      <output> stable(all) </output>
      <output> G </output>
      <output> H </output>
      <output> S </output>
      <output> Cp </output>
      <output> V </output>
      <output> dVdT </output>
      <output> dVdp </output>
      <output> G(all) </output>
      <output> H(all) </output>
      <output> S(all) </output>
      <output> V(all) </output>
    </compute>
  </ComputeOutput> 
  <ComputeOutput class=output>
    <start>
      <var name=T value=300></var>
    </start>
    <finish>
      <var name=T value=1300 operation=GE></var>
    </finish>
    <step>
      <var name=T><convert> T+100 </convert></var>
    </step>
    <compute>
      <algorithm IDREF=twosol></algorithm>
      <input value=0.5> x(YBa2Cu3O7) </input>
      <input name=T> T </input>
      <output> state(all) </output>
      <output> internal(all) </output>
      <output> stable(all) </output>
      <output> G </output>
      <output> H </output>
      <output> S </output>
      <output> Cp </output>
      <output> V </output>
      <output> dVdT </output>
      <output> dVdp </output>
      <output> G(all) </output>
      <output> H(all) </output>
      <output> S(all) </output>
      <output> V(all) </output>
    </compute>
  </ComputeOutput> 
  <ComputeOutput class=output>
    <start>
      <var name=T value=300></var>
    </start>
    <finish>
      <var name=T value=1300 operation=GE></var>
    </finish>
    <step>
      <var name=T><convert> T+100 </convert></var>
    </step>
    <compute>
      <algorithm IDREF=twosol></algorithm>
      <input once value=0.75> x(YBa2Cu3O7) </input>
      <input name=T> T </input>
      <output> state(all) </output>
      <output> internal(all) </output>
      <output> stable(all) </output>
      <output> G </output>
      <output> H </output>
      <output> S </output>
      <output> Cp </output>
      <output> V </output>
      <output> dVdT </output>
      <output> dVdp </output>
      <output> G(all) </output>
      <output> H(all) </output>
      <output> S(all) </output>
      <output> V(all) </output>
    </compute>
  </ComputeOutput> 
  <ComputeOutput class=output>
    <start>
      <var name=T value=300></var>
    </start>
    <finish>
      <var name=T value=1300 operation=GE></var>
    </finish>
    <step>
      <var name=T><convert> T+100 </convert></var>
    </step>
    <compute>
      <algorithm IDREF=twosol></algorithm>
      <input value=1> x(YBa2Cu3O7) </input>
      <input name=T> T </input>
      <output> state(all) </output>
      <output> internal(all) </output>
      <output> stable(all) </output>
      <output> G </output>
      <output> H </output>
      <output> S </output>
      <output> Cp </output>
      <output> V </output>
      <output> dVdT </output>
      <output> dVdp </output>
      <output> G(all) </output>
      <output> H(all) </output>
      <output> S(all) </output>
      <output> V(all) </output>
    </compute>
  </ComputeOutput> 
</OutputFile>

