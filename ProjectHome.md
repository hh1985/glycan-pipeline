# About HS-SEQ #
HS-SEQ is the first _de novo_ sequencing algorithm designed specifically for identifying heparan sulfate (HS) isomer using high resolution tandem mass spectrometry.  It is able to provide product-grade identification performance within seconds or minutes.

### Principle ###
HS-SEQ deals with ambiguous HS tandem MS data.  Due to the repeating disaccharide units of HS and lability of sulfate groups, peaks in the tandem MS spectra are easily to assigned to multiple fragment structures, leading to ambiguous evidence for sulfation localization. HS-SEQ constructs a graph based on the most confident peak interpretations, and expands the graph by integrating the currently most confident peak interpretations (Bayesian inference). The graph can be converted into a modification localization profile -- "modification distribution", which allows the prediction of sulfation/acetylation sites. The computation time only depends on the number of peak interpretations, and is independent of the total candidate space.

### Performance Demonstration ###
| **Test data** | **Modification distribution (dp15 for demo)** |
|:--------------|:----------------------------------------------|
| ![https://glycan-pipeline.googlecode.com/svn/trunk/GAG/image/hs_for_test.jpg](https://glycan-pipeline.googlecode.com/svn/trunk/GAG/image/hs_for_test.jpg) | ![https://glycan-pipeline.googlecode.com/svn/trunk/GAG/image/dp15_avg.jpg](https://glycan-pipeline.googlecode.com/svn/trunk/GAG/image/dp15_avg.jpg) |

### Downloads ###
Program(v1.0.1):
https://glycan-pipeline.googlecode.com/svn/trunk/GAG/bin/hsseq_v1.0.1.zip

### Usage ###
Step 0: Deisotoping/deconvolution (generating monoisotopic peak list, for Solarix only, to be uploaded) <br>
>> simplefinder.exe -s spectrum_file -p 245.000 -z 4- -o output_file -e 1 <br>

Step 1: Peak interpretation (Optional) <br>
>> librarymatch.exe -c composition_file -m monoisotopic_peak_list -o output_file -e 2 <br>
<ul><li>type "librarymatch.exe -h" for help. <br>
</li><li>composition_file example can be found from directory "structure"<br>
</li><li>monoisotopic_peak_list example can be found from directory "data" in the format: m/z intensity z</li></ul>

step 2: Sequencing <br>
>> hsseq -c composition_file -m monoisotopic_peak_list -o output_file -e 2 <br>
<ul><li>type "hsseq.exe -h" for help. <br></li></ul>

<h3>Example</h3>
>hsseq.exe -c structure\arixtra.gl -m data\arixtra_4_NETD_mo<br>
no_list.txt <br>

<a href='https://glycan-pipeline.googlecode.com/svn/trunk/GAG/image/hs-seq-output.JPG'>https://glycan-pipeline.googlecode.com/svn/trunk/GAG/image/hs-seq-output.JPG</a>

<h3>Contact</h3>
The program was developed by Han Hu(hh1985@bu.edu) <br>
PI: Prof. Joseph Zaia (jzaia@bu.edu) and Prof. Yu (Brandon) Xia (brandon.xia@mcgill.ca)<br>
<br>
<h3>Reference</h3>
Hu, H. <i>et al.</i> A Computational Framework for Heparan Sulfate Sequencing Using High-resolution Tandem Mass Spectra. <i>Mol Cell Proteomics</i> mcp.M114.039560 (2014). doi:10.1074/mcp.M114.039560 <a href='http://www.mcponline.org/content/early/2014/06/12/mcp.M114.039560.abstract'>[Link</a>]