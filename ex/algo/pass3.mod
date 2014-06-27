<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<OutputFile ext=pass3 format=file>
  <ComputeOutput>
    <start>
      <var name=x1 value=1></var>
      <var name=x2 value=10></var>
    </start>
    <finish>
      <var name=x1 value=10 operation=GE></var>
    </finish>
    <step>
      <var name=x1><convert>x1+1</convert></var>
      <var name=x2><convert>x2+10</convert></var>
    </step>
    <compute>
      <PassThrough class=algorithm id=pass></PassThrough>
      <input name=x1> x1 </input>
      <input name=x2> x2 </input>
      <output> x1 </output>
      <output> x2 </output>
      <convert name=x1>x1</convert>
      <convert name=x2>x2</convert>
      <convert name=f>exp(x1)+sin(x2)-log(x1/x2)</convert>
    </compute>
  </ComputeOutput>
</OutputFile>
