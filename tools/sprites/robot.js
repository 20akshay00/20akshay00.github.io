const {read}=require("./ase.js");const {encode}=require("./png.js");
for(const [file,tag] of [["Misc/robot_anim_new","robot"],["Void/Void","void"]]){
const a=read("C:/Users/saksh/Pictures/PixelArt/"+file+".aseprite");
a.layers.forEach((l,i)=>{const frames=[...new Set(a.cels.filter(c=>c.layer===i).map(c=>c.frame))].sort((x,y)=>x-y);console.log(tag,l.name,"frames",frames.join(","));});
// one horizontal strip per layer: every frame that layer has a cel in, on the full canvas
a.layers.forEach((l,i)=>{const cs=a.cels.filter(c=>c.layer===i).sort((x,y)=>x.frame-y.frame);if(!cs.length)return;const W=a.W,H=a.H,n=cs.length,img=Buffer.alloc(W*n*H*4);
 cs.forEach((c,k)=>{for(let y=0;y<c.h;y++)for(let x=0;x<c.w;x++){const X=c.x+x,Y=c.y+y;if(X<0||Y<0||X>=W||Y>=H)continue;c.rgba.copy(img,4*(Y*W*n+k*W+X),4*(y*c.w+x),4*(y*c.w+x)+4);}});
 encode(`sprites2/${tag}_${l.name.toLowerCase()}.png`,W*n,H,img);});}
