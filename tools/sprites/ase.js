// Minimal Aseprite reader: layers and frame-0 cels (RGBA, grayscale, indexed).
const fs=require("fs"),zlib=require("zlib");
exports.read=(file)=>{
 const b=fs.readFileSync(file);let p=0;const u8=()=>b[p++],u16=()=>{const v=b.readUInt16LE(p);p+=2;return v},s16=()=>{const v=b.readInt16LE(p);p+=2;return v},u32=()=>{const v=b.readUInt32LE(p);p+=4;return v};
 const str=()=>{const n=u16();const s=b.toString("utf8",p,p+n);p+=n;return s};
 u32();if(u16()!==0xA5E0)throw new Error("not aseprite");const frames=u16(),W=u16(),H=u16(),depth=u16();u32();u16();u32();u32();const transIdx=u8();p+=3;u16();p=128;
 const layers=[],cels=[];let palette=[];
 for(let f=0;f<frames;f++){const start=p,fbytes=u32();u16();const oldN=u16();u16();p+=2;const newN=u32();const n=newN||oldN;
  for(let c=0;c<n;c++){const cs=p,size=u32(),type=u16();
   if(type===0x2004){const flags=u16(),lt=u16(),level=u16();u16();u16();const blend=u16(),op=u8();p+=3;layers.push({name:str(),flags,type:lt,level,opacity:op,visible:!!(flags&1)});}
   else if(type===0x2019){u32();const first=u32(),last=u32();p+=8;for(let i=first;i<=last;i++){const fl=u16();palette[i]=[u8(),u8(),u8(),u8()];if(fl&1)str();}}
   else if(type===0x0004&&!palette.length){const packets=u16();let idx=0;for(let k=0;k<packets;k++){idx+=u8();let cnt=u8()||256;for(let i=0;i<cnt;i++)palette[idx++]=[u8(),u8(),u8(),255];}}
   else if(type===0x2005){const layer=u16(),x=s16(),y=s16(),op=u8(),ct=u16();s16();p+=5;if(ct===1){const lf=u16();const src=cels.find(c=>c.frame===lf&&c.layer===layer);if(src)cels.push({...src,frame:f});}
    if(ct===0||ct===2){const w=u16(),h=u16();let px=ct===0?b.subarray(p,cs+size):zlib.inflateSync(b.subarray(p,cs+size));cels.push({frame:f,layer,x,y,op,w,h,px});}}
   p=cs+size;}
  p=start+fbytes;}
 // cel pixels to RGBA
 const bpp=depth/8;
 for(const c of cels){const out=Buffer.alloc(c.w*c.h*4);for(let i=0;i<c.w*c.h;i++){
   if(depth===32)c.px.copy(out,4*i,4*i,4*i+4);
   else if(depth===16){out.fill(c.px[2*i],4*i,4*i+3);out[4*i+3]=c.px[2*i+1];}
   else{const k=c.px[i];if(k===transIdx){out[4*i+3]=0;}else{const q=palette[k]||[0,0,0,255];out[4*i]=q[0];out[4*i+1]=q[1];out[4*i+2]=q[2];out[4*i+3]=q[3];}}}
  c.rgba=out;}
 return {W,H,depth,frames,layers,cels};};
