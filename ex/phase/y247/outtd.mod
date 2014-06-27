<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<PhaseProperty class=algorithm id=Y247
  DependentMoleFraction=first>
  <phase class=phase IDREF=Y247></phase>
</PhaseProperty> 

<OutputFile ext=td format=file>
  <ComputeOutput class=output>
    <start>
      <var name=T value=300></var>
    </start>
    <finish>
      <var name=T value=1300 operation=GE></var>
    </finish>
    <step>
      <var name=T>
        <convert>
          T+100
        </convert>
      </var>
    </step>
    <compute>
      <algorithm class=algorithm IDREF=Y247></algorithm>
      <input value=0> x(Y2Ba4Cu7O15) </input>
      <input name=T> T </input>
      <output> state(all) </output>
      <output> G </output>
      <output> H </output>
      <output> S </output>
      <output> Cp </output>
      <output> G(all) </output>
      <output> H(all) </output>
      <output> S(all) </output>
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
      <var name=T>
        <convert>
          T+100
        </convert>
      </var>
    </step>
    <compute>
      <algorithm class=algorithm IDREF=Y247></algorithm>
      <input value=0.25> x(Y2Ba4Cu7O15) </input>
      <input name=T> T </input>
      <output> state(all) </output>
      <output> G </output>
      <output> H </output>
      <output> S </output>
      <output> Cp </output>
      <output> G(all) </output>
      <output> H(all) </output>
      <output> S(all) </output>
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
      <var name=T>
        <convert>
          T+100
        </convert>
      </var>
    </step>
    <compute>
      <algorithm class=algorithm IDREF=Y247></algorithm>
      <input value=0.5> x(Y2Ba4Cu7O15) </input>
      <input name=T> T </input>
      <output> state(all) </output>
      <output> G </output>
      <output> H </output>
      <output> S </output>
      <output> Cp </output>
      <output> G(all) </output>
      <output> H(all) </output>
      <output> S(all) </output>
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
      <var name=T>
        <convert>
          T+100
        </convert>
      </var>
    </step>
    <compute>
      <algorithm class=algorithm IDREF=Y247></algorithm>
      <input value=0.75> x(Y2Ba4Cu7O15) </input>
      <input name=T> T </input>
      <output> state(all) </output>
      <output> G </output>
      <output> H </output>
      <output> S </output>
      <output> Cp </output>
      <output> G(all) </output>
      <output> H(all) </output>
      <output> S(all) </output>
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
      <var name=T>
        <convert>
          T+100
        </convert>
      </var>
    </step>
    <compute>
      <algorithm class=algorithm IDREF=Y247></algorithm>
      <input value=1> x(Y2Ba4Cu7O15) </input>
      <input name=T> T </input>
      <output> state(all) </output>
      <output> G </output>
      <output> H </output>
      <output> S </output>
      <output> Cp </output>
      <output> G(all) </output>
      <output> H(all) </output>
      <output> S(all) </output>
    </compute>
  </ComputeOutput> 
</OutputFile>
