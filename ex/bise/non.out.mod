<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<OutputFile ext=non>
  <ComputeOutput class=output>
    <compute>
      <algorithm class=algorithm IDREF=Tm_s10></algorithm>
      <output> state(all) </output>
      <output> fmin </output>
      <output> fmin(all) </output>
    </compute>
  </ComputeOutput> 
  <ComputeOutput class=output>
    <compute>
      <algorithm class=algorithm IDREF=Tm_s23></algorithm>
      <output> state(all) </output>
      <output> fmin </output>
      <output> fmin(all) </output>
    </compute>
  </ComputeOutput> 
  <ComputeOutput class=output>
    <compute>
      <algorithm class=algorithm IDREF=Tm_s01></algorithm>
      <output> state(all) </output>
      <output> fmin </output>
      <output> fmin(all) </output>
    </compute>
  </ComputeOutput> 
  <ComputeOutput class=output>
    <compute>
      <algorithm class=algorithm IDREF=s10_L_s32></algorithm>
      <output> state(all) </output>
      <output> fmin </output>
      <output> fmin(all) </output>
    </compute>
  </ComputeOutput> 
  <ComputeOutput class=output>
    <compute>
      <algorithm class=algorithm IDREF=s32_L_s11></algorithm>
      <output> state(all) </output>
      <output> fmin </output>
      <output> fmin(all) </output>
    </compute>
  </ComputeOutput> 
  <ComputeOutput class=output>
    <compute>
      <algorithm class=algorithm IDREF=s11_L_s23></algorithm>
      <output> state(all) </output>
      <output> fmin </output>
      <output> fmin(all) </output>
    </compute>
  </ComputeOutput> 
    <ComputeOutput class=output>
    <compute>
      <algorithm class=algorithm IDREF=s23_L1_L2></algorithm>
      <output> state(all) </output>
      <output> fmin </output>
      <output> fmin(all) </output>
    </compute>
  </ComputeOutput> 
  <ComputeOutput class=output>
    <compute>
      <algorithm class=algorithm IDREF=s23_L_s01></algorithm>
      <output> state(all) </output>
      <output> fmin </output>
      <output> fmin(all) </output>
    </compute>
  </ComputeOutput> 
</OutputFile>
