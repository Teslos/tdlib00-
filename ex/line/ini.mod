<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<residual id=line
   yname=y
   xname=x>
  <compute>
    <PassThrough class=algorithm id=pass></PassThrough>
    <input name=x> x </input>
    <output> x </output>
    <convert>
      <coef id=a unknown=1>0</coef> + 
      <coef id=b unknown=1>0</coef>*x
    </convert>
  </compute>
</residual>

<OutputFile ext=xy format=gnuplot>
  <SeriesOutput Residuals=0>
    series_1
  </SeriesOutput>
  <SeriesOutput Residuals=0>
    series_2
  </SeriesOutput>
  <SeriesOutput Residuals=0>
    series_3
  </SeriesOutput>
  <SeriesOutput Residuals=0>
    series_4
  </SeriesOutput>
  <SeriesOutput Residuals=0>
    series_5
  </SeriesOutput>
  <SeriesOutput Residuals=0>
    series_6
  </SeriesOutput>
  <ResidualOutput InName=x OutName=y>
    <start>
      <var name=x value=0></var>
    </start>
    <finish>
      <var name=x value=100 operation=GE></var>
    </finish>
    <step>
      <var name=x><convert> x + 1 </convert></var>
    </step>
    line
  </ResidualOutput>
</OutputFile>

