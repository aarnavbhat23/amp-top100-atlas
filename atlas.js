'use strict';
const $ = id => document.getElementById(id);
const confidence = v => v>=90?'#408cff':v>=70?'#65cbf3':v>=50?'#ffdb13':'#ff7d45';
const fixed = (v,n=2) => Number(v).toFixed(n);
let rows, selected, map, molecule, models=[], pdbs=[];
function listing(id, pairs) {
  $(id).replaceChildren();
  for(const [key,value] of pairs){const dt=document.createElement('dt'),dd=document.createElement('dd');dt.textContent=key;dd.textContent=value;$(id).append(dt,dd);}
}
function styleDetail(residue){
  const kind=$('representation').value, style={};
  if(kind!=='stick')style.cartoon={colorfunc:a=>confidence(a.b)};
  if(kind!=='cartoon')style.stick={radius:.16,colorfunc:a=>confidence(a.b)};
  molecule.setStyle({},style);
  if(residue)molecule.addStyle({resi:residue},{stick:{color:'#ffffff',radius:.23}});
  molecule.render();
}
function drawPAE(r){
  const a=r.structure.pae,c=$('pae'),ctx=c.getContext('2d');ctx.clearRect(0,0,320,320);
  if(!a){$('pae-info').textContent='No PAE matrix available.';return;}
  const n=a.length,s=280/n;
  for(let y=0;y<n;y++)for(let x=0;x<n;x++){const t=Math.min(a[y][x]/30,1);ctx.fillStyle=`rgb(${Math.round(245*t)},${Math.round(70+175*t)},${Math.round(55+190*t)})`;ctx.fillRect(30+x*s,10+y*s,s+1,s+1);}
  ctx.fillStyle='#aab6c4';ctx.font='12px sans-serif';ctx.fillText('1',30,307);ctx.fillText(String(n),295,307);ctx.fillText('Residue index →',105,320);
  $('pae-info').textContent='PAE in Å: dark green = 0; white = 30 or higher. Rows and columns follow sequence order. Lower values indicate less predicted relative-position error.';
}
function select(index){
  selected=(index+rows.length)%rows.length;const r=rows[selected],s=r.structure;
  $('pick').value=selected;history.replaceState(null,'',`#top=${r.rank}`);
  $('peptide-title').textContent=`Top ${String(r.rank).padStart(3,'0')}`;
  molecule.removeAllModels();molecule.addModel(pdbs[selected],'pdb');styleDetail();molecule.zoomTo();molecule.render();
  $('structure-status').textContent=`${s.model} · sequence verified · predicted monomer`;
  $('plddt').textContent=fixed(s.mean_ca_plddt,1);$('ptm').textContent=fixed(s.metrics.ptm,3);$('length').textContent=r.length;
  $('sequence').replaceChildren();
  [...r.sequence].forEach((aa,i)=>{const b=document.createElement('button');b.textContent=aa;b.style.color=confidence(s.residue_plddt[i]);b.title=`${aa}${i+1}: pLDDT ${fixed(s.residue_plddt[i],1)}`;b.onclick=()=>{styleDetail(i+1);$('residue-info').textContent=b.title;};$('sequence').append(b);});
  $('residue-info').textContent='Select a residue to highlight its atoms.';
  for(const key of ['pdb','cif'])$(key).href=s[key];$('metrics').href=s.metrics_file;
  const c=r.calculated;
  listing('properties',[['Average molecular weight',`${fixed(c.molecular_weight_da,1)} Da`],['Charge at pH 7',fixed(c.charge_pH7)],['K + R − D − E',r.charge_proxy_KR_DE],['Predicted pI (algorithm estimate)',fixed(c.pI)],['Hydrophobic fraction (AVILMFWY)',fixed(r.hydrophobic_fraction_AVILMFWY,3)],['GRAVY (Kyte–Doolittle)',fixed(c.gravy,3)],['Aromatic fraction (FYW)',fixed(c.aromatic_fraction,3)],['Cysteines',r.cysteines],['Longest identical-residue run',r.longest_homopolymer],['Composition',Object.entries(c.counts).filter(([,n])=>n).map(([a,n])=>`${a}:${n}`).join(' ')]]);
  listing('reference',[['Maximum reference ratio',fixed(r.max_reference_ratio,4)],['Nearest reference',r.nearest_reference_id],['Reference sequence',r.nearest_reference_sequence],['Maximum within-top100 ratio',fixed(r.max_internal_ratio,4)],['Exact DAMP membership',r.damp_exact_sources.join(', ')||'None in the audited DAMP snapshot'],['Library membership',r.present_in_library==='True'?'Confirmed':'Not confirmed']]);
  $('warnings').textContent=Number(r.max_reference_ratio)>=.8?'At the 0.8000 boundary: the audited validator rejects >0.8, not ≥0.8.': 'Passes audited sequence constraints. This is not certification of future challenge rules or experimental activity.';
  if(Number(r.cysteines))$('warnings').textContent+=' One cysteine; no intramolecular disulfide is implied.';
  $('job').textContent=`Job ${s.job_id}. Original CIF retained; PDB uses residue pLDDT in its B-factor field.`;drawPAE(r);
}
function recolour(){
  const key=$('colour').value;
  const values=rows.map(r=>key==='confidence'?r.structure.mean_ca_plddt:key==='charge'?r.calculated.charge_pH7:key==='similarity'?+r.max_reference_ratio:+r[key]);
  const lo=Math.min(...values),hi=Math.max(...values);
  models.forEach((m,i)=>{const t=(values[i]-lo)/(hi-lo||1);const rgb=[Math.round(64+191*t),Math.round(140+17*t),Math.round(255-190*t)];const col=key==='confidence'?confidence(values[i]):'#'+rgb.map(v=>v.toString(16).padStart(2,'0')).join('');m.setStyle({},{cartoon:{color:col}});});map.render();
  $('map-legend').textContent=key==='confidence'?'Mean pLDDT: blue ≥90 · cyan 70–90 · yellow 50–70 · orange <50':`${fixed(lo)} (blue) → ${fixed(hi)} (orange)`;
}
async function main(){
  const response=await fetch('atlas_complete.json');if(!response.ok)throw Error('Data could not be loaded');const data=await response.json();rows=data.peptides;
  pdbs=await Promise.all(rows.map(async r=>{const res=await fetch(r.structure.pdb);if(!res.ok)throw Error(`Missing structure ${r.rank}`);return res.text();}));
  map=$3Dmol.createViewer($('map'),{backgroundColor:'#080c13'});molecule=$3Dmol.createViewer($('molecule'),{backgroundColor:'#101720'});
  const center=[0,1,2].map(k=>rows.reduce((sum,r)=>sum+r.umap[k],0)/rows.length);
  rows.forEach((r,i)=>{
    const lines=pdbs[i].split('\n'), atoms=lines.filter(l=>l.startsWith('ATOM  '));
    const centroid=[30,38,46].map(k=>atoms.reduce((sum,l)=>sum+Number(l.slice(k,k+8)),0)/atoms.length);
    const translated=lines.map(l=>{if(!l.startsWith('ATOM  ')&&!l.startsWith('HETATM'))return l;const xyz=[30,38,46].map((k,j)=>(Number(l.slice(k,k+8))-centroid[j]+(r.umap[j]-center[j])*65).toFixed(3).padStart(8));return l.slice(0,30)+xyz.join('')+l.slice(54);}).join('\n');
    const m=map.addModel(translated,'pdb');models.push(m);
    m.setClickable({},true,()=>select(i));m.setHoverable({},true,()=>{$('hover').textContent=`Top ${r.rank} · ${r.sequence} · mean pLDDT ${fixed(r.structure.mean_ca_plddt,1)}`;},()=>{$('hover').textContent='';});
    const option=document.createElement('option');option.value=i;option.textContent=`${String(r.rank).padStart(3,'0')} · ${r.sequence}`;$('pick').append(option);
  });
  recolour();map.zoomTo();map.render();select(Math.max(0,rows.findIndex(r=>+r.rank===Number(location.hash.split('=')[1]||1))));
  $('count').textContent=`${data.structure_summary.ready}/100 predicted structures · all sequences verified`;
  $('pick').onchange=()=>select(+$('pick').value);$('previous').onclick=()=>select(selected-1);$('next').onclick=()=>select(selected+1);
  $('search').oninput=()=>{const q=$('search').value.trim().toUpperCase();if(!q){$('search-status').textContent='';return;}const matches=rows.map((r,i)=>({r,i})).filter(({r})=>/^\d+$/.test(q)?+r.rank===+q:r.sequence.includes(q));$('search-status').textContent=`${matches.length} match(es)`;if(matches.length)select(matches[0].i);};
  $('zoom-in').onclick=()=>{map.zoom(1.3);map.render();};$('zoom-out').onclick=()=>{map.zoom(.77);map.render();};$('reset-map').onclick=()=>{map.zoomTo();map.render();};$('focus').onclick=()=>{map.zoomTo({model:models[selected].getID()});map.render();};
  $('colour').onchange=recolour;$('representation').onchange=()=>styleDetail();$('reset-molecule').onclick=()=>{molecule.zoomTo();molecule.render();};
  $('expand').onclick=()=>{document.body.classList.toggle('expanded');$('expand').textContent=document.body.classList.contains('expanded')?'Back to atlas':'Expand structure';setTimeout(()=>{map.resize();molecule.resize();molecule.zoomTo();molecule.render();},100);};
  $('copy').onclick=async()=>{try{await navigator.clipboard.writeText(rows[selected].sequence);$('copy').textContent='Copied';setTimeout(()=>$('copy').textContent='Copy sequence',1200);}catch{$('residue-info').textContent=rows[selected].sequence;}};
  $('methods').onclick=()=>$('method-dialog').showModal();$('close-methods').onclick=()=>$('method-dialog').close();$('explore').onclick=()=>$('method-dialog').close();
  const paragraphs=[
    '100 peptides from Vinay’s top100.csv, source commit f77bfdbbf2aaaa4442d7453591d8a1c6ee224588. Submission rank is retained; this app does not re-rank candidates.',
    'Boltz-2.1 single-chain monomer predictions. Every downloaded structure was checked for an exact sequence match, one chain, and finite coordinates. Native CIF and raw metrics are downloadable. Confidence is not antibacterial activity or experimental validation. pTM and pLDDT measure different aspects of confidence; neither establishes a stable solution conformation for a short peptide.',
    'ESM-2 t12 35M, final-layer residue-mean embeddings (480 dimensions). UMAP: cosine metric, 15 neighbors, min_dist 0.2, 3 dimensions, seed 42. Only these 100 peptides are embedded. Ribbon glyphs are translated to map positions, not structurally aligned; the detail viewer retains original coordinates.',
    `Layout diagnostics at k=10: trustworthiness ${fixed(data.method.trustworthiness_cosine_k10,3)}, original-neighbor recall ${fixed(data.method.mean_original_neighbor_recall_k10,3)}, seed-42/43 neighbor overlap ${fixed(data.method.mean_seed_neighbor_overlap_k10,3)}. Apparent clusters are exploratory, not evidence of functional families.`,
    'Reference audit: 39,448 sequences; 3,944,800 top100/reference comparisons using the validator’s normalized Levenshtein ratio, not alignment percent identity. No ratio >0.8; top40 equals 0.8. The audited full library contains 50,000 unique canonical sequences. This validates the supplied export, not the undocumented generation/ranking pipeline.',
    'Sequence properties use Biopython 1.88 with unmodified free termini. pI is an algorithmic estimate (boundary values near 12 should not be treated as precise). Exact DAMP membership is a snapshot lookup, not proof of novelty or training exposure.',
    'Local visualization: 3Dmol.js 2.5.3. No structures are experimental; no membrane or binding partner was modeled. All 100 API cost estimates summed to $2.50; actual billing was not verified.'
  ];for(const t of paragraphs){const p=document.createElement('p');p.textContent=t;$('method-content').append(p);}
  window.addEventListener('resize',()=>{map.resize();molecule.resize();});
}
main().catch(e=>{$('count').textContent=`Unable to load atlas: ${e.message}. Serve this folder through a local HTTP server.`;console.error(e);});
