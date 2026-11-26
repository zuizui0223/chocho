// Final layout only. Preserve source-v07 artwork and all prediction outcomes.
const fs=require("fs"),path=require("path");
const root=process.argv[2]||".";
const render=function renderFinal(f2,f4){
  function sub(s,x,y){if(s.split(x).length!==2)throw Error("Nonunique anchor "+x.slice(0,120));return s.replace(x,y)}
  for(const [x,y] of [["+0.2854","+0.29"],["+0.0954","+0.10"],["+0.1900","+0.19"],["0.046030","0.046"],["0.004487","0.0045"]])f2=sub(f2,">"+x+"</text>",">"+y+"</text>"); 
  f2=sub(f2,"Introduced ≈10.26× native median","Introduced ≈10× native median");
  let start=f2.indexOf('<rect x="20" y="1131"');let end=f2.indexOf('<rect x="555" y="1131"');
  if(start<0||end<=start)throw Error("E panel markers absent");
  let e=f2.slice(start,end);
  e=e.replaceAll("#bd673c","#2e6876");
  const oldLegend=e.indexOf('<text x="35" y="1383"');
  if(oldLegend<0)throw Error("legend anchor absent");
  const label=(x,y,s,weight=400,size=19)=>`<text x="${x}" y="${y}" font-family="Arial,Helvetica,sans-serif" font-size="${size}" fill="#4f6473" font-weight="${weight}">${s}</text>\n`;
  e=e.slice(0,oldLegend)+`<circle cx="43" cy="1383" r="5.8" fill="#2e6876"/><circle cx="158" cy="1383" r="5" fill="#fff" stroke="#2e6876" stroke-width="2"/>\n`+label(57,1389,"Base",700)+label(172,1389,"MEM10",700)+label(272,1389,"95% cluster intervals",400,17)+label(35,1414,"Species − family: p = 0.008",400,18)+label(35,1437,"With MEM10: p = 0.059",400,18);
  e=sub(e,'height="279"','height="313"');
  f2=f2.slice(0,start)+e+f2.slice(end);
  f2=sub(f2,'<rect x="555" y="1131" width="525" height="279"','<rect x="555" y="1131" width="525" height="313"');
  f2=sub(f2,'height="1425" viewBox="0 0 1100 1425"','height="1459" viewBox="0 0 1100 1459"');
  // Keep A,B exactly as in the verified six-panel source; C uses its original graphic.
  const c0=f4.indexOf('<rect x="20" y="618" width="573" height="353"');
  const d0=f4.indexOf('<rect x="607" y="618" width="573" height="353"');
  const e0=f4.indexOf('<rect x="20" y="984" width="573" height="353"');
  const f0=f4.indexOf('<rect x="607" y="984" width="573" height="353"');
  const end4=f4.indexOf('<text x="23" y="1370"');
  if(!(c0>0&&c0<d0&&d0<e0&&e0<f0&&f0<end4))throw Error("Panel boundaries");
  const txt=(x,y,s,size=19,weight=400,color="#172c40")=>`<text x="${x}" y="${y}" font-family="Arial,Helvetica,sans-serif" font-size="${size}" font-weight="${weight}" fill="${color}">${s}</text>\n`;
  let lead=f4.slice(0,c0);
  lead=sub(lead,'height="1435" viewBox="0 0 1200 1435"','height="1005" viewBox="0 0 1200 1005"');
  lead=sub(lead,'height="1435" rx="0"','height="1005" rx="0"');
  lead=sub(lead,"C–F: original independent prediction scores and adverse transfer.","C: original independent later-record prediction scores. Other prediction tests are in Figure S5.2.");
  // Re-use exact C chart rather than alter any data, shift center; increase its height in overall page by leaving bottom margin.
  let main4=lead+`<g transform="translate(303.5 0)">\n`+f4.slice(c0,d0)+`</g>\n`+txt(22,992,"A–B: mapped hosts and indexed records are not colonization evidence. C: retrospective predictions; further tests in Fig. S5.2.",17,400,"#4f6473")+"</svg>\n";
  // Preserve D/E/F source geometry, labels, scores, and whiskers; only relocate panels and relabel A/B/C.
  const labelPanel=(s,from,to)=>sub(s,">"+from+"</text>",">"+to+"</text>");
  const top=`<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="830" viewBox="0 0 1200 830" role="img" aria-labelledby="title desc">\n<title id="title">Additional butterfly community prediction and adverse transfer checks</title>\n<desc id="desc">Previous main Figure 4D–F, relettered A–C without changing numerical evidence. A: community prediction; B: unseen divisions; C: adverse South American transfer.</desc>\n<rect width="1200" height="830" fill="white"/>\n`+txt(24,36,"S5.2  |  Recorded community prediction and geographic transfer",27,700)+txt(25,64,"Relative MSE reductions; intervals are conditional Monte Carlo quantiles, not population confidence intervals.",19,400,"#4f6473");
  const supp=top+`<g transform="translate(-587 -533)">\n`+labelPanel(f4.slice(d0,e0),"D","A")+`</g>\n<g transform="translate(587 -899)">\n`+labelPanel(f4.slice(e0,f0),"E","B")+`</g>\n<g transform="translate(-283.5 -533)">\n`+labelPanel(f4.slice(f0,end4),"F","C")+`</g>\n`+txt(24,818,"Positive: improved held-out MSE; negative: worse predictions. Full South America contrasts are retained in Fig. S5.1.",18,400,"#4f6473")+"</svg>\n";
  return [f2,main4,supp];
};
const f2=fs.readFileSync(path.join(root,"figures/integrated_island_port_butterfly_v03/Fig2_ports_native_alien_floras_readable_v07.svg"),"utf8");
const f4=fs.readFileSync(path.join(root,"figures/integrated_island_port_butterfly_v03/Fig4_recorded_butterfly_isolation_monarch_v07.svg"),"utf8");
const outputs=render(f2,f4), names=["figures/integrated_island_port_butterfly_v03/Fig2_ports_native_alien_floras_readable_v08.svg","figures/integrated_island_port_butterfly_v03/Fig4_recorded_butterfly_isolation_monarch_v08.svg","figures/integrated_island_port_butterfly_v03/FigS5_2_butterfly_community_prediction_transfer_v08.svg"];
for(let i=0;i<names.length;i++){const p=path.join(root,names[i]);if(process.argv.includes("--check")){if(fs.readFileSync(p,"utf8")!==outputs[i])throw Error("Final figure mismatch "+p);}else fs.writeFileSync(p,outputs[i]);}
console.log("V08_MAIN_AND_SI_RENDERINGS_MATCH_SOURCE");
