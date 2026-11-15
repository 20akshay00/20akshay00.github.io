// minimal PNG decode (8-bit RGBA / RGB / palette) and RGBA encode
const zlib=require("zlib"),fs=require("fs");
function crcTable(){const t=new Int32Array(256);for(let n=0;n<256;n++){let c=n;for(let k=0;k<8;k++)c=c&1?0xedb88320^(c>>>1):c>>>1;t[n]=c;}return t;}
const T=crcTable();const crc=(b)=>{let c=-1;for(const x of b)c=T[(c^x)&255]^(c>>>8);return (c^-1)>>>0;};
exports.decode=(file)=>{const b=fs.readFileSync(file);let p=8,W,H,bd,ct,idat=[],plte=null,trns=null;
 while(p<b.length){const len=b.readUInt32BE(p),type=b.toString("ascii",p+4,p+8),d=b.slice(p+8,p+8+len);
  if(type==="IHDR"){W=d.readUInt32BE(0);H=d.readUInt32BE(4);bd=d[8];ct=d[9];}else if(type==="IDAT")idat.push(d);else if(type==="PLTE")plte=d;else if(type==="tRNS")trns=d;p+=12+len;}
 if(bd!==8)throw new Error("bit depth "+bd);
 const ch={2:3,6:4,3:1,0:1,4:2}[ct],bpp=ch,stride=W*bpp,raw=zlib.inflateSync(Buffer.concat(idat));
 const px=Buffer.alloc(W*H*bpp);let prev=Buffer.alloc(stride);
 for(let y=0;y<H;y++){const f=raw[y*(stride+1)],line=raw.subarray(y*(stride+1)+1,(y+1)*(stride+1)),out=Buffer.alloc(stride);
  for(let i=0;i<stride;i++){const a=i>=bpp?out[i-bpp]:0,up=prev[i],c=i>=bpp?prev[i-bpp]:0;let v=line[i];
   if(f===1)v+=a;else if(f===2)v+=up;else if(f===3)v+=(a+up)>>1;else if(f===4){const q=a+up-c,pa=Math.abs(q-a),pb=Math.abs(q-up),pc=Math.abs(q-c);v+=(pa<=pb&&pa<=pc)?a:(pb<=pc?up:c);}out[i]=v&255;}
  out.copy(px,y*stride);prev=out;}
 const rgba=Buffer.alloc(W*H*4);
 for(let i=0;i<W*H;i++){if(ct===6)px.copy(rgba,4*i,4*i,4*i+4);else if(ct===2){rgba[4*i]=px[3*i];rgba[4*i+1]=px[3*i+1];rgba[4*i+2]=px[3*i+2];rgba[4*i+3]=255;}
  else if(ct===3){const k=px[i];rgba[4*i]=plte[3*k];rgba[4*i+1]=plte[3*k+1];rgba[4*i+2]=plte[3*k+2];rgba[4*i+3]=trns&&k<trns.length?trns[k]:255;}
  else if(ct===0){rgba.fill(px[i],4*i,4*i+3);rgba[4*i+3]=255;}else{rgba.fill(px[2*i],4*i,4*i+3);rgba[4*i+3]=px[2*i+1];}}
 return {W,H,rgba};};
exports.encode=(file,W,H,rgba)=>{const raw=Buffer.alloc(H*(W*4+1));for(let y=0;y<H;y++){raw[y*(W*4+1)]=0;rgba.copy(raw,y*(W*4+1)+1,y*W*4,(y+1)*W*4);}
 const chunk=(t,d)=>{const l=Buffer.alloc(4);l.writeUInt32BE(d.length);const td=Buffer.concat([Buffer.from(t,"ascii"),d]);const c=Buffer.alloc(4);c.writeUInt32BE(crc(td));return Buffer.concat([l,td,c]);};
 const ih=Buffer.alloc(13);ih.writeUInt32BE(W,0);ih.writeUInt32BE(H,4);ih[8]=8;ih[9]=6;
 fs.writeFileSync(file,Buffer.concat([Buffer.from([137,80,78,71,13,10,26,10]),chunk("IHDR",ih),chunk("IDAT",zlib.deflateSync(raw,{level:9})),chunk("IEND",Buffer.alloc(0))]));};
