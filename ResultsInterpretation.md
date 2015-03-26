# Introduction #

The result of HS-SEQ is called "Modification distribution"
https://glycan-pipeline.googlecode.com/svn/trunk/GAG/image/mod_dist.JPG

# Explanation #
  * Residue: type of monosaccharide residue.<br>
<ul><li>Position: Residue ID (from non-reducing end to reducing end) - Site ID. For example, GlcA 1-2 means 2-<i>O</i> modification on the 2nd (not the 1st) residue which is a GlcA.  GlcN 6-6 means 6-<i>O</i> modification on the 7th residue which is a GlcN<br>
</li><li>Modification: acetylation or sulfation.<br>
</li><li>Intensity: likelihood for modification.</li></ul>

<h1>Exception</h1>
Sulfate loss or ambiguous peak interpretation may lead to weird modification distribution, but this can be tell easily. <br>
Example: The following distribution suggests GlcN(ID: 0) is highly likely to have one sulfate group on it (prob = 0.326364 x 3 = 0.979), while GlcA (ID: 1) is unlikely to have one (prob = 0.0294). However, since no cross-ring information is used on GlcN, all candidate sites (2-<i>N</i>, 3-<i>O</i>, 6-<i>O</i>) have equal chance to be sulfated <br>
GlcN	0-2	SO3	0.326364 <br>
GlcN	0-3	SO3	0.326364 <br>
GlcN	0-6	SO3	0.326364 <br>
GlcA	1-2	SO3	0.0293662 <br>

<h1>Notice</h1>
If a cross-ring cleavage causes significant sulfate loss, the probibility of sulfated sites on the cleaved residue may be wrong, but the total number of sulfate groups on the residue is very likely to be correct. That's one of the advantages of HS-SEQ: error is localized.