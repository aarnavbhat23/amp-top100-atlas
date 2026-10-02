export const aminoAcids='ACDEFGHIKLMNPQRSTVWY';
const kd={A:1.8,C:2.5,D:-3.5,E:-3.5,F:2.8,G:-.4,H:-3.2,I:4.5,K:-3.9,L:3.8,M:1.9,N:-3.5,P:-1.6,Q:-3.5,R:-4.5,S:-.8,T:-.7,V:4.2,W:-.9,Y:-1.3};
const mw={A:89.0932,C:121.1582,D:133.1027,E:147.1293,F:165.1891,G:75.0666,H:155.1546,I:131.1729,K:146.1876,L:131.1729,M:149.2113,N:132.1179,P:115.1305,Q:146.1445,R:174.201,S:105.0926,T:119.1192,V:117.1463,W:204.2252,Y:181.1885};
export function features(sequence){
 const counts=Object.fromEntries([...aminoAcids].map(a=>[a,0]));for(const a of sequence)counts[a]++;
 const n=sequence.length, sum=(letters)=>[...letters].reduce((s,a)=>s+counts[a],0);
 const entropy=Object.values(counts).reduce((s,c)=>c?s-c/n*Math.log2(c/n):s,0)/Math.log2(20);
 return {length:n,charge_proxy:sum('KR')-sum('DE'),hydrophobic_fraction:sum('AVILMFWY')/n,gravy:[...sequence].reduce((s,a)=>s+kd[a],0)/n,molecular_weight:[...sequence].reduce((s,a)=>s+mw[a],0)-18.0153*(n-1),aromatic_fraction:sum('FYW')/n,cysteines:counts.C,entropy,counts};
}
export function delimited(text,delimiter){
 const rows=[];let row=[],field='',quoted=false;
 for(let i=0;i<text.length;i++){const c=text[i];if(c==='"'){if(quoted&&text[i+1]==='"'){field+='"';i++;}else if(quoted){quoted=false;}else if(!field){quoted=true;}else throw Error('Unexpected quote in an unquoted field.');}
 else if(c===delimiter&&!quoted){row.push(field);field='';}else if((c==='\n'||c==='\r')&&!quoted){if(c==='\r'&&text[i+1]==='\n')i++;row.push(field);if(row.some(v=>v.trim()))rows.push(row);row=[];field='';}else field+=c;}
 if(quoted)throw Error('Unclosed quoted field.');row.push(field);if(row.some(v=>v.trim()))rows.push(row);return rows;
}
export function parseFile(text,name){
 text=text.replace(/^\uFEFF/,'');
 if(text.trimStart().startsWith('>')||/\.(fa|faa|fasta|fna)$/i.test(name)){
 const records=[];let current=null;for(const line of text.split(/\r?\n/)){if(!line.trim())continue;if(line.startsWith('>')){current={id:line.slice(1).trim(),raw:''};records.push(current);}else{if(!current)throw Error('FASTA sequence appears before its header.');current.raw+=line.trim();}}
 return {kind:'fasta',records};}
 const table=delimited(text,/\.tsv$/i.test(name)?'\t':',');if(table.length<2)throw Error('CSV/TSV needs a header and at least one record.');
 const headers=table[0];if(table.slice(1).some(r=>r.length!==headers.length))throw Error('Inconsistent column counts. Check the delimiter and quoted fields.');
 return {kind:'table',headers,rows:table.slice(1)};
}
export function analyze(records){
 const seen=new Set();return records.map((r,i)=>{
 const sequence=r.raw.replace(/\s/g,'').toUpperCase(),invalid=[...new Set([...sequence].filter(a=>!aminoAcids.includes(a)))].join('');
 const valid=sequence.length>0&&!invalid,duplicate=valid&&seen.has(sequence);if(valid)seen.add(sequence);
 const issues=[];if(!sequence)issues.push('Empty sequence');if(invalid)issues.push(`Invalid letters: ${invalid}`);if(valid&&(sequence.length<8||sequence.length>50))issues.push('Outside challenge length 8–50');if(duplicate)issues.push('Duplicate sequence');
 return {row:i+1,id:r.id||`row-${i+1}`,sequence,valid,duplicate,normalized:r.raw!==sequence,issues,features:valid?features(sequence):null};});
}
