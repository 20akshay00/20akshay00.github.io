const {read}=require("./ase.js");const {encode}=require("./png.js");const root="C:/Users/saksh/Pictures/PixelArt/";
const files=["Scenes/WizardFloat","Scenes/MechPhi","Scenes/escape","Scenes/trapped","Scenes/heroic-Recovered","Scenes/light","Scenes/darkness","Scenes/Halves","Scenes/Bot1","Scenes/Char","Scenes/Wizard1","Scenes/Phi@i-Recovered","Scenes/Phi@i","Scenes/cloud","Scenes/choices","Dragon/duality","Misc/Robot","Void/Void","Void/face","Void/lamp","SecretSanta2023/player","SecretSanta2023/attractor","GWJ59/master","GWJ59/disciple","LD55/card","idgbjam/terminal","RadarDefender/base","GMTK2023/target"];
const skip=/^(background1?|bg|layer 1\(hidden\))$/i;
for(const f of files){const a=read(root+f+".aseprite");const W=a.W,H=a.H,img=Buffer.alloc(W*H*4);
 a.layers.forEach((l,i)=>{if(!l.visible||l.type===1||skip.test(l.name))return;
  for(const c of a.cels.filter(c=>c.layer===i)){const op=(c.op/255)*(l.opacity/255);
   for(let y=0;y<c.h;y++)for(let x=0;x<c.w;x++){const X=c.x+x,Y=c.y+y;if(X<0||Y<0||X>=W||Y>=H)continue;const s=4*(y*c.w+x),t=4*(Y*W+X),sa=c.rgba[s+3]/255*op;if(!sa)continue;const da=img[t+3]/255,oa=sa+da*(1-sa);
    for(let k=0;k<3;k++)img[t+k]=Math.round((c.rgba[s+k]*sa+img[t+k]*da*(1-sa))/oa);img[t+3]=Math.round(oa*255);}}});
 let x0=W,y0=H,x1=-1,y1=-1;for(let y=0;y<H;y++)for(let x=0;x<W;x++)if(img[4*(y*W+x)+3]>8){x0=Math.min(x0,x);x1=Math.max(x1,x);y0=Math.min(y0,y);y1=Math.max(y1,y);}
 if(x1<0){console.log(f,"empty");continue;}const w=x1-x0+1,h=y1-y0+1,out=Buffer.alloc(w*h*4);for(let y=0;y<h;y++)img.copy(out,y*w*4,4*((y0+y)*W+x0),4*((y0+y)*W+x0+w));
 const name=f.split("/").pop().replace(/[@ ]/g,"_");encode(`sprites2/${name}.png`,w,h,out);console.log(name,w,h,"layers:",a.layers.filter((l,i)=>l.visible&&!skip.test(l.name)&&l.type!==1).map(l=>l.name).join(","));}
