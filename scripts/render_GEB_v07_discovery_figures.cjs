// Deterministic editorial rendering of verified display inputs. No model selection.
const fs=require("fs"),path=require("path");
const root=process.argv[2]||".";
const folder=path.join(root,"figures/integrated_island_port_butterfly_v03/");
const render=function render(d,o,src){
 const blue="#2e6876", orange="#bd673c",ink="#172c40",sub="#4f6473",bg="#f5f8fa",grid="#dce4e8",white="#fff";
 const esc=s=>String(s).replace(/&/g,"&amp;").replace(/</g,"&lt;").replace(/>/g,"&gt;");
 const num=v=>Number(v.toFixed(3)),pct=(v,dig=1)=>(100*v).toFixed(dig)+"%",sgn=(v,dig=2)=>(v>=0?"+":"−")+Math.abs(v).toFixed(dig);
 const t=(x,y,s,size=19,color=ink,bold=false,anchor="start",italic=false)=>`<text x="${num(x)}" y="${num(y)}" font-family="Arial,Helvetica,sans-serif" font-size="${size}" fill="${color}" font-weight="${bold?700:400}" text-anchor="${anchor}"${italic?' font-style="italic"':""}>${esc(s)}</text>\n`;
 const rect=(x,y,w,h,fill=bg,rx=8)=>`<rect x="${num(x)}" y="${num(y)}" width="${num(w)}" height="${num(h)}" rx="${rx}" fill="${fill}"/>\n`;
 const line=(x,y,x2,y2,color=grid,sw=1,dash="")=>`<line x1="${num(x)}" y1="${num(y)}" x2="${num(x2)}" y2="${num(y2)}" stroke="${color}" stroke-width="${sw}"${dash?` stroke-dasharray="${dash}"`:""}/>\n`;
 const circ=(x,y,r,c,stroke="none",opacity=1)=>`<circle cx="${num(x)}" cy="${num(y)}" r="${num(r)}" fill="${c}" stroke="${stroke}" stroke-width="1.8" opacity="${opacity}"/>\n`;
 const header=(x,y,w,h,letter,title,subtitle="")=>rect(x,y,w,h)+t(x+16,y+32,letter,24,blue,true)+t(x+53,y+32,title,21,ink,true)+(subtitle?t(x+16,y+58,subtitle,17,sub):"");
 const append=(svg,h0,h1,extra,desc)=>{
   let s=svg.replace(new RegExp(`height="${h0}" viewBox="0 0 (\\d+) ${h0}"`),(_,w)=>`height="${h1}" viewBox="0 0 ${w} ${h1}"`);
   if(s===svg)throw Error("root dimensions did not change");
   s=s.replace(/<desc id="desc">([\s\S]*?)<\/desc>/,(_,v)=>`<desc id="desc">${v} ${esc(desc)}</desc>`);
   return s.replace(/<\/svg>\s*$/,extra+"</svg>\n");
 };
 let f1=o[0].replace('height="284.4pt" viewBox="0 0 720 284.4"','height="316.4pt" viewBox="0 0 720 316.4"');
 f1=f1.replace(/<\/svg>\s*$/,rect(28.8,285,677,27,"#e9f1f1",3)+t(40,303,`Introduced-only share of current opportunities:  ${pct(d.island_intro_share)} islands  |  ${pct(d.mainland_intro_share)} mainlands`,12.7,ink,true)+"</svg>\n");
 // Figure 2: keep all original panels, add hierarchy and trait comparison.
 let e=header(20,1131,525,279,"E","Host identity matters","126 islands · standardized distance slopes");
 const ex=v=>190+(v+.35)/.52*332;
 for(const v of [-.3,-.2,-.1,0,.1]){e+=line(ex(v),1200,ex(v),1339,v===0?sub:grid,v===0?1.8:1,v===0?"4 4":"");e+=t(ex(v),1359,v===0?"0":sgn(v,1),17,sub,false,"middle")}
 d.hierarchy.forEach((a,i)=>{const y=1217+45*i;e+=t(35,y+6,a[0],19,ink,true);e+=line(ex(a[2]),y-6,ex(a[3]),y-6,blue,3)+circ(ex(a[1]),y-6,5.8,blue);e+=line(ex(a[5]),y+9,ex(a[6]),y+9,orange,2.2)+circ(ex(a[4]),y+9,5,white,orange)});
 e+=t(35,1383,"Base",17,blue,true)+t(102,1383,"○ MEM10",17,orange,true)+t(224,1383,"Nine-Level1 95% intervals",17,sub);
 e+=t(35,1402,"Species − family: p = 0.008; MEM10 + cluster p = 0.059",15.8,sub);
 let ff=header(555,1131,525,279,"F","Seawater-dispersal morphology","22 regions · native minus introduced share");
 let fx=v=>590+v/6*460;
 for(let v of [0,2,4,6]){ff+=line(fx(v),1206,fx(v),1290,grid);ff+=t(fx(v),1311,String(v),17,sub,false,"middle")}
 ff+=line(fx(d.seawater.descriptive_interval[0]),1250,fx(d.seawater.descriptive_interval[1]),1250,blue,4)+circ(fx(d.seawater.difference_pp),1250,7,blue);
 ff+=t(fx(d.seawater.difference_pp),1231,"+2.95 pp",22,blue,true,"middle")+t(817,1333,"Native − introduced (percentage points)",17,sub,false,"middle");
 ff+=t(574,1360,"17 positive · 4 negative · 1 tie",20,ink,true)+t(574,1385,"Descriptive resampling interval; biased trait coverage.",17,sub);
 ff+=t(574,1403,"Morphology is not observed current-driven dispersal.",16,sub);
 let f2=append(o[1],1120,1425,e+ff,"E: base and MEM10 taxonomic rank distance slopes with clustered intervals. F: native-minus-introduced seawater dispersal morphology.");
 // Figure 3: preserve original A-C, add explicit fraction of initially empty opportunities.
 let p3=header(23,878,1173,215,"D","More of the gaps filled for broadly supplied butterflies","Pooled island gaps filled, ordered by global native-resource extent (not diet breadth)");
 d.gap_groups.forEach((a,i)=>{let x=57+i*287,v=a[1]/a[2];p3+=t(x+116,969,["Narrowest","Second","Third","Broadest"][i]+` · n=${a[0]}`,18,ink,true,"middle")+rect(x,980,232,25,"#dfe6e9",4)+rect(x,980,232*v/.5,25,i===0?blue:orange,4)+t(x+116,1032,pct(v),25,i===0?blue:orange,true,"middle")+t(x+116,1056,`${a[1].toLocaleString("en-US")} / ${a[2].toLocaleString("en-US")} gaps`,17,sub,false,"middle")});
 p3+=t(42,1082,"Gap filling is pooled added / initially empty units (not gain relative to native coverage).",17,sub);
 let f3=append(o[2],870,1106,p3,"D: pooled resource gaps filled in the four groups by global native host geography.");
 // Figure 4: new A-B first, old four numerical contrasts re-rendered C-F.
 const W=1200,H=1435;
 let f4=`<svg xmlns="http://www.w3.org/2000/svg" width="${W}" height="${H}" viewBox="0 0 ${W} ${H}" role="img" aria-labelledby="title desc">\n<title id="title">Independent butterfly records and mapped host resources</title>\n<desc id="desc">A: raw recorded share points and source-refit adjusted recorded and potential shares for the same 84 regions. B: monarch introduced-only host geography and records, with a named endemic feeding example. C–F: original independent prediction scores and adverse transfer.</desc>\n`+rect(0,0,W,H,white,0)+t(24,36,"4  |  From mapped host opportunities to recorded butterflies",27,ink,true)+t(25,65,"Separate observational and predictive outcomes; no claim of colonization or demographic gain.",19,sub);
 // A.
 const x=20,y=85,w=573,h=520;
 f4+=header(x,y,w,h,"A","Observed records and resource isolation","224 shared butterflies · 84 record-active regions");
 const ax=85,ay=180,aw=472,ah=249;
 const xp=v=>ax+Math.log1p(v)/Math.log1p(6000)*aw,yp=v=>ay+ah*(1-v);
 for(const v of [0,.25,.5,.75,1])f4+=line(ax,yp(v),ax+aw,yp(v))+t(ax-10,yp(v)+6,pct(v,0),18,sub,false,"end");
 for(const v of [0,10,100,1000,6000])f4+=line(xp(v),ay,xp(v),ay+ah)+t(xp(v),ay+ah+24,String(v),18,sub,false,"middle");
 f4+=line(ax,ay+ah,ax+aw,ay+ah,sub,1.5)+t(ax+aw/2,ay+ah+48,"Mainland distance (km; log1p scale)",18,sub,false,"middle");
 d.observations.forEach(([c,dist,n,i])=>{f4+=circ(xp(dist),yp(i/(n+i)),2.5+Math.sqrt(n+i)*.65,orange,"#fff",.55)});
 for(const k of ["potential","observed"]){const [mu,b]=d.fits[k];let points=[];for(let i=0;i<=60;i++){let l=Math.log1p(d.fits.min_km)+(Math.log1p(d.fits.max_km)-Math.log1p(d.fits.min_km))*i/60;let q=1/(1+Math.exp(-(mu+b*(l-d.fits.mean_logdistance)/d.fits.sd_logdistance)));points.push(`${num(xp(Math.expm1(l)))},${num(yp(q))}`)}f4+=`<polyline points="${points.join(" ")}" fill="none" stroke="${k==="potential"?blue:orange}" stroke-width="4" ${k==="potential"?'stroke-dasharray="8 5"':""}/>\n`}
 f4+=line(113,154,149,154,blue,4,"7 4")+t(157,160,"Mapped potential fit",17,blue,true)+line(328,154,364,154,orange,4)+t(372,160,"Recorded fit",17,orange,true);
 f4+=t(35,502,"Nearest quartile",18,sub,true)+t(35,532,"11.4%  (70/612)",25,orange,true)+t(35,555,"32/32 regions represented",17,sub);
 f4+=t(305,502,"Farthest quartile",18,sub,true)+t(305,532,"61.4%  (27/44)",25,orange,true)+t(305,555,"12/32 regions represented",17,sub);
 f4+=t(35,580,"Adjusted recorded β = +0.44; clustered p = 0.040",18,ink,true)+t(35,599,"Effort-restricted results uncertain; dots sized by record-unit count.",16,sub);
 // B.
 const bx=607,by=85,bw=573,bh=520;
 f4+=header(bx,by,bw,bh,"B","Monarch: introduced-only host regions","42 botanical regions; 20 with GBIF records");
 const mx=625,my=168,mw=537,mh=mw/360*145;
 f4+=rect(mx,my,mw,mh,"#eaf1f3",0);
 f4+=`<g transform="translate(${mx} ${my}) scale(${mw/360})"><path d="${d.basemap.path}" fill="#cdd7da" stroke="#a7b8bd" stroke-width=".45"/></g>\n`;
 for(const lon of [-120,-60,0,60,120])f4+=line(mx+(lon+180)*mw/360,my,mx+(lon+180)*mw/360,my+mh,"#c7d6d9",.8,"3 3");
 for(const lat of [-30,0,30,60])f4+=line(mx,my+(85-lat)*mw/360,mx+mw,my+(85-lat)*mw/360,"#c7d6d9",.8,"3 3");
 for(const [c,lon,lat,record] of [...d.monarch].sort((a,b)=>a[3]-b[3])){let X=mx+(lon+180)*mw/360,Y=my+(85-lat)*mw/360;f4+=`<g><title>${esc(c)}: ${record?"recorded":"not indexed in this source"}</title>${circ(X,Y,4.8,record?orange:white,record?white:sub,.94)}</g>`}
 const haw=d.monarch.find(x=>x[0]==="HAW"), hx=mx+(haw[1]+180)*mw/360,hy=my+(85-haw[2])*mw/360;
 f4+=circ(hx,hy,9,"none",ink)+line(hx+7,hy-6,hx+36,hy-25,ink,1.4)+t(hx+39,hy-28,"Hawaii",18,ink,true);
 f4+=circ(637,410,5,orange,white)+t(650,416,"Recorded (20)",18,orange,true)+circ(829,410,5,white,sub)+t(842,416,"Not indexed (22)",18,sub);
 f4+=t(625,442,"18 of 20 recorded in multiple years.",19,ink,true)+t(625,466,"1875 account: monarch appearance in Hawaii",18,sub)+t(625,488,"followed milkweed introduction (Scudder).",18,sub);
 f4+=line(625,505,1162,505,grid)+t(625,531,"A separate endemic feeding example",19,blue,true)+t(625,556,"Pieris cheiranthi larvae on introduced Tropaeolum majus;",17,ink)+t(625,579,"Tenerife and La Palma; native Crambe also documented.",17,sub)+t(625,599,"Region markers, not record coordinates or transport routes.",16,sub);
 // Lower four panels, exact receipt numbers, interval types kept separate.
 const p=src.figures.fig4;
 function bars(x,y,w,h,L,title,subtitle,items,lo,hi,ticks,note1,note2){
   let s=header(x,y,w,h,L,title,subtitle);
   let left=x+198,right=x+w-64,top=y+84,row=54,bottom=y+239;
   let xx=v=>left+(v-lo)/(hi-lo)*(right-left);
   for(const v of ticks){s+=line(xx(v),top-8,xx(v),bottom,v===0?sub:grid,v===0?1.7:1);s+=t(xx(v),bottom+23,`${v}%`,17,sub,false,"middle")}
   items.forEach((it,i)=>{let yy=top+22+i*row;let v=it[1];s+=t(x+15,yy+6,it[0],18,ink);s+=rect(Math.min(xx(0),xx(v)),yy-10,Math.max(1,Math.abs(xx(v)-xx(0))),20,v<0?blue:orange,2);if(it[2]){s+=line(xx(it[2][0]),yy,xx(it[2][1]),yy,ink,2)+line(xx(it[2][0]),yy-6,xx(it[2][0]),yy+6,ink,2)+line(xx(it[2][1]),yy-6,xx(it[2][1]),yy+6,ink,2)}s+=t(x+w-10,yy+6,sgn(v,it[3]??2)+"%",17,v<0?blue:orange,true,"end")});
   s+=t(x+15,y+h-41,note1,17,sub)+t(x+15,y+h-18,note2,16,sub);return s;
 }
 const aa=p.panel_A_155_independently_indexed_butterflies.mse_reduction_percent,bb=p.panel_B_224_butterfly_source,cc=p.panel_C_eight_unseen_worlds,dd=p.panel_D_south_america;
 f4+=bars(20,618,573,353,"C","Prediction of later records","155 butterflies · 208,482 test cases",[["Temporal",aa.time_held,null,2],["Spatial block",aa.spatial_held,null,2],["≥100 km distant",aa.over_100km,null,2]],0,5,[0,2.5,5],"Relative MSE reduction (host information added).","≥100 km is a subset, not an independent fit.");
 f4+=bars(607,618,573,353,"D","Recorded community prediction","224 butterflies · 169 regions · 14,196 pairs",[["Original",bb.original_H_after_O_MSE_improvement_pct,null,3],["Dataset-days",bb.dataset_days_MSE_improvement_pct_median,bb.dataset_days_MC95,3],["Record counts",bb.record_count_MSE_improvement_pct_median,bb.record_count_MC95,3]],0,.8,[0,.4,.8],"Known-host flora added after other alien flora.","Whiskers: conditional MC quantiles (29 draws).");
 f4+=bars(20,984,573,353,"E","Transfer to unseen divisions","Eight Level1 divisions · world-equal loss",[["Original",cc.original_H_after_O_world_equal_pct,null,3],["Dataset-days",cc.dataset_day_median,cc.dataset_day_MC95,3],["Record counts",cc.record_count_median,cc.record_count_MC95,3]],0,1.8,[0,.9,1.8],"Positive average; not universal across divisions.","Whiskers: conditional MC quantiles (29 draws).");
 f4+=bars(607,984,573,353,"F","South America: adverse transfer","2,980 cross-world pairs · 20-region holdout",[["Dataset-days",dd.dataset_day_cross_median,dd.dataset_day_cross_MC95,2],["Record counts",dd.record_count_cross_median,dd.record_count_cross_MC95,2]],-8.5,0,[-8,-4,0],"Within-SA: −0.496% / −0.126% (190 pairs).","Full four-stratum comparison remains in SI S5.1.");
 f4+=t(23,1370,"A: back-transformed equal-region fitted log odds at mean covariates; raw dots and quartile totals use recorded units.",18,sub);
 f4+=t(23,1394,"B: simplified GSHHG coastline; mapped current hosts and 2001–2024 indices are not dated colonization evidence.",18,sub);
 f4+=t(23,1418,"C–F: within-model relative MSE reductions; negative values indicate poorer predictions. Do not pool populations.",18,sub)+"</svg>\n";
 return [f1,f2,f3,f4];
};
const d=JSON.parse(fs.readFileSync(path.join(folder,"FIGURE_V07_DISPLAY_SOURCE.json"),"utf8"));
const originals=["Fig1_resource_geography_scales_v03.svg","Fig2_ports_native_alien_floras_readable_v06.svg","Fig3_global_host_contributors_at_ports_restructured_v04.svg"].map(n=>fs.readFileSync(path.join(folder,n),"utf8"));
const s=JSON.parse(fs.readFileSync(path.join(folder,"FIGURE_V04_SOURCE_NUMBERS.json"),"utf8"));
const figures=render(d,originals,s);const names=["Fig1_resource_geography_scales_v07.svg","Fig2_ports_native_alien_floras_readable_v07.svg","Fig3_global_host_contributors_at_ports_v07.svg","Fig4_recorded_butterfly_isolation_monarch_v07.svg"];
for(let i=0;i<names.length;i++){let p=path.join(folder,names[i]);if(process.argv.includes("--check")){if(fs.readFileSync(p,"utf8")!==figures[i])throw Error("Figure source mismatch "+p);}else fs.writeFileSync(p,figures[i]);}
console.log("FOUR_V07_RENDERINGS_MATCH_SOURCE");
