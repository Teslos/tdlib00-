<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<OutputFile ext=non format=file>
  <ComputeOutput class=output>
    <compute>
      <algorithm class=algorithm IDREF=TmBa></algorithm>
      <output> state(all) </output>
      <output> fmin </output>
      <output> fmin(all) </output>
    </compute>
  </ComputeOutput> 
  <ComputeOutput class=output>
    <compute>
      <algorithm class=algorithm IDREF=Ba_L_BaCu></algorithm>
      <output> state(all) </output>
      <output> fmin </output>
      <output> fmin(all) </output>
    </compute>
  </ComputeOutput> 
  <ComputeOutput class=output>
    <compute>
      <algorithm class=algorithm IDREF=BaCu_L_BaCu13></algorithm>
      <output> state(all) </output>
      <output> fmin </output>
      <output> fmin(all) </output>
    </compute>
  </ComputeOutput> 
  <ComputeOutput class=output>
    <compute>
      <algorithm class=algorithm IDREF=BaCu13_L_Cu></algorithm>
      <output> state(all) </output>
      <output> fmin </output>
      <output> fmin(all) </output>
    </compute>
  </ComputeOutput> 
  <ComputeOutput class=output>
    <compute>
      <algorithm class=algorithm IDREF=TmCu></algorithm>
      <output> state(all) </output>
      <output> fmin </output>
      <output> fmin(all) </output>
    </compute>
  </ComputeOutput> 
</OutputFile>

<OutputFile ext=mis format=file>
  <SpinodalOutput phase=L direction=up></SpinodalOutput>
  <SpinodalOutput phase=L direction=down></SpinodalOutput>
</OutputFile>
