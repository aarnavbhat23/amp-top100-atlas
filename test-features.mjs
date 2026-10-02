import assert from 'node:assert/strict';
import fs from 'node:fs';
import {parseFile,analyze,features} from './sequence-features.mjs';
assert.deepEqual(parseFile('>one\nac de\nFG\n>two\nKKKKKKKK\n','a.fasta').records,[{id:'one',raw:'ac deFG'},{id:'two',raw:'KKKKKKKK'}]);
assert.equal(parseFile('id,sequence\r\n"a,b",ACDEFGHI\r\n','a.csv').rows[0][0],'a,b');
assert.equal(parseFile('id\tsequence\n1\tACDEFGHI','a.tsv').rows[0][1],'ACDEFGHI');
assert.throws(()=>parseFile('id,seq\n"bad,KK','a.csv'));
assert.throws(()=>parseFile('id,seq\na,KK,extra','a.csv'));
const a=analyze([{id:'a',raw:' acdefghi '},{id:'b',raw:'ACDEFGHI'},{id:'c',raw:'ABC*'},{id:'d',raw:''}]);
assert(a[0].valid&&a[0].normalized);assert(a[1].duplicate);assert(!a[2].valid&&!a[3].valid);
assert.equal(features('KKDDEERR').charge_proxy,0);assert.equal(features('AAAAAAAA').entropy,0);
const data=JSON.parse(fs.readFileSync('atlas_complete.json')).peptides;
for(const r of data){const f=features(r.sequence);assert.equal(f.length,+r.length);assert.equal(f.charge_proxy,+r.charge_proxy_KR_DE);assert(Math.abs(f.gravy-r.calculated.gravy)<1e-10);assert(Math.abs(f.molecular_weight-r.calculated.molecular_weight_da)<1e-6);assert(Math.abs(f.hydrophobic_fraction-Number(r.hydrophobic_fraction_AVILMFWY))<1e-10);}
console.log('Parsing, validation, duplicate tests passed; descriptors agree with all 100 existing audited peptides.');
