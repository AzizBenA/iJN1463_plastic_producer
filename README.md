# IJN1436 Plastic Producer

# iJN1463 Updates — Plastics / C2 Platform (Project Log)

<p align="center">
  <b>Repository note:</b> This README tracks incremental edits to <code>iJN1463</code> (updated variants) to support
  C2 platform constraints and PET/PU/PBAT-derived monomer assimilation (EG, TA, AA, BDO), with curated reactions,
  bounds edits, and balancing status.
</p>

<hr/>

## 📌 Quick summary

<ul>
  <li><b>C2 platform integration:</b> key bounds constrained (e.g., <code>ACS</code>) and new acetate phosphorylation reaction (<code>ACPH</code>).</li>
  <li><b>PET branch:</b> terephthalate metabolites + dioxygenase/dehydrogenase reactions added.</li>
  <li><b>PU/PBAT branch:</b> ethylene glycol, BDO, 4-hydroxybutyrate path added; reactions balanced.</li>
  <li><b>Yield “fix” constraints:</b> selected reverse-TCA and threonine-to-acetyl-CoA routes must be blocked (plasmid-free strain context).</li>
  <li><b>Open item:</b> plastics “polymer” metabolites were added but still need proper <b>charge annotation</b>.</li>
</ul>

<hr/>

## 🧱 Model variants and change log

### <code>iJN1463_updated_1</code> — C2 platform integration
<ul>
  <li><b>C2 platform present</b></li>
  <li><code>ACS</code> bounds set to zero</li>
  <li><code>ALDD2y</code>, <code>ALCD2x</code> bounds set to zero</li>
  <li><code>ACPH</code> created: acetate phosphorylation</li>
  <li><code>3OXCOAT</code> bounds set to zero</li>
  <li>New reactions: <code>ADPCOAR</code>, <code>ADPCOAH</code></li>
  <li>New metabolites: adipic acid <code>adpac_c</code>, adipoyl-CoA <code>adpcoa_c</code></li>
</ul>

<hr/>

### <code>250211_iJN1463_updated_2</code> — PET monomers: terephthalate path
<ul>
  <li>Added metabolite: terephthalic acid <code>terepa_c</code></li>
  <li>Added metabolite: 1,2-dihydroxy-1,2-dihydroterephthalate <code>12di12ditere_c</code></li>
  <li>Added reaction: <code>TEREDEOXY</code> (Terephthalate 1,2-dioxygenase)</li>
  <li>Added reaction: <code>TEREDEHYD</code> (Terephthalate 1,2-dehydrogenase)</li>
</ul>

<hr/>

### <code>250324_iJN1463_updated_3</code> — PU/PBAT monomers: EG & BDO module + 4HB branch
<ul>
  <li>Added metabolite: Ethylene glycol <code>etheglycol_c</code></li>
  <li>Added reaction: <code>ETHYGLYCHYD</code> (Ethylene glycol dehydrogenase)</li>

  <li>Added metabolite: Butane-1,4-diol <code>14btdl_c</code></li>
  <li>Added metabolite: 4-hydroxybutyraldehyde <code>4hbutald_c</code></li>
  <li>Added reaction: <code>14BUTADEH</code> (1,4-butanediol dehydrogenase)</li>

  <li>Added metabolite: 4-hydroxybutyrate <code>h4but_c</code></li>
  <li>Added reaction: <code>4HYDBUDEH</code> (4-hydroxybutyraldehyde dehydrogenase)</li>
  <li>Added reaction: <code>4HYDBUSUCDEH</code> (4-hydroxybutyrate → succinic semialdehyde)</li>

  <li>Added reaction: <code>AACOA3HYD</code> (Acetoacetyl-CoA → (R)-3-hydroxybutyryl-CoA)</li>

  <li><code>AACOAT</code>, <code>BDH</code> bounds set to zero</li>
  <li><b>Balancing status:</b> all reactions mass/charge balanced except (see next version)</li>
</ul>

<hr/>

### <code>250325_iJN1463_updated_4</code> — Balancing complete + utilization note
<ul>
  <li><b>Balancing status:</b> all reactions are balanced</li>
  <li><b>Note:</b> the model is not using <code>(R)-3-Hydroxybutanoate</code>; “you need Hao” (pending follow-up)</li>
</ul>

<hr/>

### <code>250709_iJN1463_updated_5</code> — Bounds edits / KO set for expected yield
<p>Knockouts / constraints applied to match the “right yield” behavior:</p>

```python
model.reactions.get_by_id('SUCOAS').bounds = (-1000.0, 0.0)
model.reactions.get_by_id('PPC').bounds    = (-0.0, 1000.0)
model.reactions.get_by_id('PC').bounds     = (-0.0, 0.0)
model.reactions.get_by_id('THRS').bounds   = (-0.0, 0.0)
