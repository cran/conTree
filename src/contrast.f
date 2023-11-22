c     mortran 2.0     (version of 7/04/75 mod 7/4/87 (ajc))             
      subroutine fcontrast (no,ni,x,y,y2,z,w,lx,mxt,itre,rtre,mxc,cat,ms
     *,isc)
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      integer lx(ni),ms(no,ni,2),itre(6,mxt),isc(no)                    
      double precision x(no,ni),y(no),y2(no),z(no),w(no),rtre(4,mxt),cat
     *(mxc)
      save nodes,nct,kstor,lstor                                        
      data mxtrm /10/                                                   
      call dosort(no,ni,x,w,ms,nu,isc)                                  
      lstor=0                                                           
      itre(4,1)=0                                                       
      itre(5,1)=1                                                       
      itre(6,1)=nu                                                      
      itre(2,1)=2                                                       
      itre(3,1)=3                                                       
      call andarm(no,y,y2,z,w,rtre(2,1),rtre(4,1))                      
      nct=1                                                             
      call split7(no,ni,x,y,y2,z,w,lx,ms,1,nu,itre(1,1),rtre(1,1),  cri1
     *,cri2,w1,w2,rtre(3,1),kct,cat(nct))
      if(itre(1,1) .ne. 0)goto 10021                                    
      itre(4,1)=-9999                                                   
      return                                                            
10021 continue                                                          
      if(itre(1,1) .ge. 0)goto 10041                                    
      rtre(1,1)=nct                                                     
      nct=nct+kct                                                       
      if(nct.gt.mxc) lstor=1                                            
10041 continue                                                          
      rtre(2,2)=0.0                                                     
      rtre(2,3)=rtre(2,2)                                               
      rtre(3,2)=cri1                                                    
      rtre(3,3)=cri2                                                    
      itre(4,2)=-1                                                      
      itre(4,3)=itre(4,2)                                               
      rtre(4,2)=w1                                                      
      rtre(4,3)=w2                                                      
      rtre(2,1)=sign(max(0d0,rtre(3,1)-rtre(2,1)),rtre(3,1))            
      if(lstor.ne.0) return                                             
      nodes=3                                                           
      ntrm=0                                                            
      if(rtre(2,1).gt.0.0) ntrm=1                                       
      ktrg=1                                                            
      kstor=0                                                           
      continue                                                          
10051 if(ntrm.ge.mxtrm)goto 10052                                       
      jt=itre(1,ktrg)                                                   
      st=rtre(1,ktrg)                                                   
      k5=itre(5,ktrg)                                                   
      k6=itre(6,ktrg)                                                   
      ju=0                                                              
      kp=0                                                              
      kl=0                                                              
      kr=0                                                              
      if(jt .ge. 0)goto 10071                                           
      ju=-jt                                                            
      kp=int(st+0.1)                                                    
      np=int(abs(cat(kp))+0.1)                                          
10071 continue                                                          
      do 10081 j=1,ni                                                   
      if(ni.gt.1.and.j.eq.jt)goto 10081                                 
      kl=k5-1                                                           
      kr=k6+1                                                           
      do 10091 i=k5,k6                                                  
      l=ms(i,j,1)                                                       
      if(jt .le. 0)goto 10111                                           
      if(x(l,jt) .ge. st)goto 10131                                     
      kl=kl+1                                                           
      isc(kl)=l                                                         
      goto 10141                                                        
10131 continue                                                          
      kr=kr-1                                                           
      isc(kr)=l                                                         
10141 continue                                                          
      continue                                                          
      goto 10151                                                        
10111 continue                                                          
      in=0                                                              
      do 10161 im=1,np                                                  
      if(x(l,ju).ne.cat(kp+im))goto 10161                               
      in=1                                                              
      goto 10162                                                        
10161 continue                                                          
10162 continue                                                          
      if(cat(kp) .le. 0.0)goto 10181                                    
      if(in .ne. 0)goto 10201                                           
      kl=kl+1                                                           
      isc(kl)=l                                                         
      goto 10211                                                        
10201 continue                                                          
      kr=kr-1                                                           
      isc(kr)=l                                                         
10211 continue                                                          
      continue                                                          
      goto 10221                                                        
10181 continue                                                          
      if(in .eq. 0)goto 10241                                           
      kl=kl+1                                                           
      isc(kl)=l                                                         
      goto 10251                                                        
10241 continue                                                          
      kr=kr-1                                                           
      isc(kr)=l                                                         
10251 continue                                                          
      continue                                                          
10221 continue                                                          
      continue                                                          
10151 continue                                                          
      continue                                                          
10091 continue                                                          
      continue                                                          
      do 10261 i=k5,kl                                                  
      ms(i,j,1)=isc(i)                                                  
10261 continue                                                          
      continue                                                          
      do 10271 i=kr,k6                                                  
      ms(i,j,1)=isc(k6-i+kr)                                            
10271 continue                                                          
      continue                                                          
10081 continue                                                          
      continue                                                          
      itre(5,itre(2,ktrg))=k5                                           
      itre(6,itre(2,ktrg))=kl                                           
      itre(5,itre(3,ktrg))=kr                                           
      itre(6,itre(3,ktrg))=k6                                           
      nde=itre(2,ktrg)                                                  
      cri=0d0                                                           
      call split7(no,ni,x,y,y2,z,w,lx,ms,itre(5,nde),itre(6,nde),itre(1,
     *nde),  rtre(1,nde),cri1,cri2,w1,w2,cri,kct,cat(nct))
      if(itre(1,nde) .ge. 0)goto 10291                                  
      rtre(1,nde)=nct                                                   
      nct=nct+kct                                                       
      if(nct.gt.mxc) lstor=1                                            
10291 continue                                                          
      if(nodes+2 .le. mxt)goto 10311                                    
      kstor=1                                                           
      goto 10321                                                        
10311 continue                                                          
      itre(2,nde)=nodes+1                                               
      itre(3,nde)=nodes+2                                               
      l=nodes+1                                                         
      rtre(2,l)=0.0                                                     
      rtre(2,l+1)=rtre(2,l)                                             
      rtre(3,l)=cri1                                                    
      rtre(3,l+1)=cri2                                                  
      rtre(4,l)=w1                                                      
      rtre(4,l+1)=w2                                                    
      itre(4,l)=-nde                                                    
      itre(4,l+1)=itre(4,l)                                             
      rtre(2,nde)=sign(max(0d0,cri-rtre(3,nde)),cri)                    
10321 continue                                                          
      continue                                                          
      nde=itre(3,ktrg)                                                  
      call split7(no,ni,x,y,y2,z,w,lx,ms,itre(5,nde),itre(6,nde),itre(1,
     *nde),  rtre(1,nde),cri1,cri2,w1,w2,cri,kct,cat(nct))
      if(itre(1,nde) .ge. 0)goto 10341                                  
      rtre(1,nde)=nct                                                   
      nct=nct+kct                                                       
      if(nct.gt.mxc) lstor=1                                            
10341 continue                                                          
      if(nodes+4 .le. mxt)goto 10361                                    
      kstor=1                                                           
      goto 10371                                                        
10361 continue                                                          
      itre(2,nde)=nodes+3                                               
      itre(3,nde)=nodes+4                                               
      l=nodes+3                                                         
      rtre(2,l)=0.0                                                     
      rtre(2,l+1)=rtre(2,l)                                             
      rtre(3,l)=cri1                                                    
      rtre(3,l+1)=cri2                                                  
      rtre(4,l)=w1                                                      
      rtre(4,l+1)=w2                                                    
      itre(4,l)=-nde                                                    
      itre(4,l+1)=itre(4,l)                                             
      rtre(2,nde)=sign(max(0d0,cri-rtre(3,nde)),cri)                    
      itre(4,ktrg)=-itre(4,ktrg)                                        
      rsv=rtre(2,ktrg)                                                  
      crix=0.0                                                          
      do 10381 k=1,nodes                                                
      if(itre(4,k).ge.0)goto 10381                                      
      if(abs(rtre(2,k)).le.crix)goto 10381                              
      if(itre(1,k).eq.0)goto 10381                                      
      crix=abs(rtre(2,k))                                               
      ktrg=k                                                            
10381 continue                                                          
      continue                                                          
10371 continue                                                          
      continue                                                          
      if(crix.le.0.0.or.kstor.ne.0.or.lstor.ne.0) return                
      nodes=nodes+4                                                     
      if(rsv.gt.0.0) ntrm=ntrm+1                                        
      goto 10051                                                        
10052 continue                                                          
      return                                                            
      entry set_trm(irg)                                                
      mxtrm=max0(irg,2)                                                 
      return                                                            
      entry get_stor(irg1,irg2)                                         
      irg1=nodes                                                        
      irg2=nct-1                                                        
      return                                                            
      entry get_err(irg1,irg2)                                          
      irg1=kstor                                                        
      irg2=lstor                                                        
      return                                                            
      end                                                               
      subroutine ans (x,itre,rtre,cat,yh)                               
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      integer itre(6,*)                                                 
      double precision x(*),rtre(4,*),cat(*)                            
      k=1                                                               
      continue                                                          
10391 if(itre(4,k).lt.0)goto 10392                                      
      if(itre(1,k) .le. 0)goto 10411                                    
      if(x(itre(1,k)) .ge. rtre(1,k))goto 10431                         
      k=itre(2,k)                                                       
      goto 10441                                                        
10431 continue                                                          
      k=itre(3,k)                                                       
10441 continue                                                          
      continue                                                          
      goto 10451                                                        
10411 continue                                                          
      j=-itre(1,k)                                                      
      kp=int(rtre(1,k)+0.1)                                             
      in=0                                                              
      np=int(abs(cat(kp))+0.1)                                          
      do 10461 i=1,np                                                   
      if(x(j).ne.cat(kp+i))goto 10461                                   
      in=1                                                              
      goto 10462                                                        
10461 continue                                                          
10462 continue                                                          
      if(cat(kp) .le. 0.0)goto 10481                                    
      if(in .ne. 0)goto 10501                                           
      k=itre(2,k)                                                       
      goto 10511                                                        
10501 continue                                                          
      k=itre(3,k)                                                       
10511 continue                                                          
      continue                                                          
      goto 10521                                                        
10481 continue                                                          
      if(in .eq. 0)goto 10541                                           
      k=itre(2,k)                                                       
      goto 10551                                                        
10541 continue                                                          
      k=itre(3,k)                                                       
10551 continue                                                          
      continue                                                          
10521 continue                                                          
      continue                                                          
10451 continue                                                          
      continue                                                          
      goto 10391                                                        
10392 continue                                                          
      yh=rtre(3,k)                                                      
      return                                                            
      end                                                               
      subroutine dosort(no,ni,x,w,ms,nu,isc)                            
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      integer ms(no,ni,*),isc(no)                                       
      double precision x(no,ni),w(no)                                   
      data new /1/                                                      
      isc = isc + 0                                                     
      if(new .eq. 0)goto 10571                                          
      do 10581 j=1,ni                                                   
      do 10591 i=1,no                                                   
      ms(i,j,1)=i                                                       
10591 continue                                                          
      continue                                                          
      call psort8(x(1,j),ms(1,j,1),1,no)                                
10581 continue                                                          
      continue                                                          
      do 10601 j=1,ni                                                   
      do 10611 i=1,no                                                   
      ms(i,j,2)=ms(i,j,1)                                               
10611 continue                                                          
      continue                                                          
      nu=0                                                              
      do 10621 i=1,no                                                   
      if(w(ms(i,j,2)).le.0.0)goto 10621                                 
      nu=nu+1                                                           
      ms(nu,j,1)=ms(i,j,2)                                              
10621 continue                                                          
      continue                                                          
10601 continue                                                          
      continue                                                          
      return                                                            
10571 continue                                                          
      do 10631 j=1,ni                                                   
      nu=0                                                              
      do 10641 i=1,no                                                   
      if(w(ms(i,j,2)).le.0.0)goto 10641                                 
      nu=nu+1                                                           
      ms(nu,j,1)=ms(i,j,2)                                              
10641 continue                                                          
      continue                                                          
10631 continue                                                          
      continue                                                          
      return                                                            
      entry set_new(irg)                                                
      new=irg                                                           
      return                                                            
      end                                                               
      subroutine split7 (no,ni,x,y,y2,z,w,lx,m,m1,m2,  jt,sp,cri1,cri2,w
     *1,w2,crm,kct,cat)
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      parameter(maxcat=1000, big=9.9e35)                                
      integer lx(ni),m(no,ni)                                           
      double precision x(no,ni),y(no),y2(no),z(no),w(no),cat(*),tcat(max
     *cat)
      data xmiss,ntn,pwr /9.0e35,500,2/                                 
      crx=0.0                                                           
      crm=0.0                                                           
      jt=0                                                              
      if(m2-m1+1.lt.2*ntn) return                                       
      do 10651 j=1,ni                                                   
      if(x(m(m1,j),j).ge.x(m(m2,j),j).or.lx(j).eq.0)goto 10651          
      if(lx(j) .ne. 1)goto 10671                                        
      call eav(x(1,j),y,y2,z,w,m(1,j),m1,m2,ntn,pwr,xmiss,tp,  cril,criu
     *,wl,wu,cri)
      if(cri.eq.-xmiss)goto 10651                                       
      if(abs(cri) .lt. crm)goto 10691                                   
      crm=abs(cri)                                                      
      cri1=cril                                                         
      cri2=criu                                                         
      w1=wl                                                             
      w2=wu                                                             
      jt=j                                                              
      sp=tp                                                             
      crx=cri                                                           
10691 continue                                                          
      goto 10701                                                        
10671 continue                                                          
      call ceav(x(1,j),y,y2,z,w,m(1,j),m1,m2,ntn,pwr,lct,  tcat,cril,cri
     *u,wl,wu,cri)
      if(cri .lt. crm)goto 10721                                        
      crm=cri                                                           
      cri1=cril                                                         
      cri2=criu                                                         
      w1=wl                                                             
      w2=wu                                                             
      jt=-j                                                             
      kct=lct                                                           
      do 10731 k=1,kct                                                  
      cat(k)=tcat(k)                                                    
10731 continue                                                          
      continue                                                          
10721 continue                                                          
10701 continue                                                          
      continue                                                          
10651 continue                                                          
      continue                                                          
      if(jt.gt.0) crm=crx                                               
      return                                                            
      entry set_miss(arg)                                               
      xmiss=arg                                                         
      return                                                            
      entry set_ntn(irg)                                                
      ntn=max0(irg,3)                                                   
      return                                                            
      entry set_pwr(arg)                                                
      pwr=arg                                                           
      return                                                            
      end                                                               
      subroutine eav (x,y,y2,z,w,m,m1,m2,ntn,pwr,xmiss,tp,  cri1s,cri2s,
     *w1s,w2s,cri)
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      integer m(*)                                                      
      double precision x(*),y(*),y2(*),z(*),w(*)                        
      data nint,icri /4,1/                                              
      xl=x(m(m1))                                                       
      if(xl .lt. xmiss)goto 10751                                       
      cri=-xmiss                                                        
      return                                                            
10751 continue                                                          
      nu=m2                                                             
      k=m(nu)                                                           
      continue                                                          
10761 if(x(k).lt.xmiss)goto 10762                                       
      nu=nu-1                                                           
      k=m(nu)                                                           
      goto 10761                                                        
10762 continue                                                          
      crx=0.0                                                           
      crimx=crx                                                         
      if(nu .ge. m2)goto 10781                                          
      call andarm(nu-m1+1,y(m(m1:nu)),y2(m(m1:nu)),  z(m(m1:nu)),w(m(m1:
     *nu)),cri1xs,w1x)
      call andarm(m2-nu,y(m((nu+1):m2)),y2(m((nu+1):m2)),  z(m((nu+1):m2
     *)),w(m((nu+1):m2)),cri2xs,w2x)
      crx=max(cri1xs,cri2xs)                                            
      crimx=(crx**pwr)*float(nu-m1+1)*float(m2-nu)/float(m2-m1+1)**2    
10781 continue                                                          
      rq=(nu-m1+1)                                                      
      mq=int(rq/nint)                                                   
      kq=1                                                              
      crim=-xmiss                                                       
      cris=-xmiss                                                       
      do 10791 i=nu,m1+1,-1                                             
      k=m(i)                                                            
      tq=0.5*(x(k)+x(m(i-1)))                                           
      if(tq.le.x(m(i-1)))goto 10791                                     
      if(tq.ge.x(k))goto 10791                                          
      if(i-m1.lt.ntn.or.nu-i+1.lt.ntn)goto 10791                        
      if(i .ge. nu-kq*mq+1)goto 10811                                   
      kq=kq+1                                                           
      call andarm(i-m1,y(m(m1:(i-1))),y2(m(m1:(i-1))),  z(m(m1:(i-1))),w
     *(m(m1:(i-1))),cri1,w1)
      call andarm(nu-i+1,y(m(i:nu)),y2(m(i:nu)),  z(m(i:nu)),w(m(i:nu)),
     *cri2,w2)
      if(icri .ne. 1)goto 10831                                         
      cri=max(cri1,cri2)                                                
      goto 10841                                                        
10831 continue                                                          
      cri=abs(cri1-cri2)                                                
10841 continue                                                          
      continue                                                          
      crit=(cri**pwr)*float(i-m1)*float(nu-i+1)/float(nu-m1+1)**2       
      if(crit.lt.crim)goto 10791                                        
      crim=crit                                                         
      tp=tq                                                             
      cris=cri                                                          
      w1s=w1                                                            
      w2s=w2                                                            
      cri1s=cri1                                                        
      cri2s=cri2                                                        
10811 continue                                                          
10791 continue                                                          
      continue                                                          
      cri=cris                                                          
      if(x(m(m2)).lt.xmiss) return                                      
      tp=xmiss                                                          
      cri1s=cri1xs                                                      
      cri2s=cri2xs                                                      
      w1s=w1x                                                           
      w2s=w2x                                                           
      if(crx .le. cri)goto 10861                                        
      cri=crx                                                           
      return                                                            
10861 continue                                                          
      cri=-cri                                                          
      return                                                            
      entry set_qint(irg)                                               
      nint=irg                                                          
      return                                                            
      entry set_cri(irg)                                                
      icri=irg                                                          
      return                                                            
      end                                                               
      subroutine ceav (x,y,y2,z,w,m,m1,m2,ntn,pwr,kct,cat,cri1,cri2,w1,w
     *2,cri)
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      parameter(maxcat=1000)                                            
      integer m(*),mt(maxcat)                                           
      double precision x(*),y(*),y2(*),z(*),w(*),cat(*),v(maxcat,3)     
      fntn=ntn                                                          
      l=0                                                               
      i1=m1                                                             
      do 10871 i=m1,m2-1                                                
      k=m(i)                                                            
      if(x(m(i+1)).le.x(k))goto 10871                                   
      l=l+1                                                             
      v(l,1)=x(k)                                                       
      i2=i-1                                                            
      call andarm(i2-i1+1,y(m(i1:i2)),y2(m(i1:i2)),  z(m(i1:i2)),w(m(i1:
     *i2)),v(l,2),v(l,3))
      i1=i                                                              
10871 continue                                                          
      continue                                                          
      k=m(m2)                                                           
      l=l+1                                                             
      v(l,1)=x(k)                                                       
      call andarm(m2-i1+1,y(m(i1:m2)),y2(m(i1:m2)),z(m(i1:m2)),w(m(i1:m2
     *)),  v(l,2),v(l,3))
      do 10881 i=1,l                                                    
      mt(i)=i                                                           
10881 continue                                                          
      continue                                                          
      call psort8(v(1:l,2),mt,1,l)                                      
      do 10891 j=1,l                                                    
      v(j,2)=v(j,2)*v(j,3)                                              
10891 continue                                                          
      continue                                                          
      sl=0.0                                                            
      wl=sl                                                             
      cri=wl                                                            
      scri=cri                                                          
      sr=sum(v(1:l,2))                                                  
      wr=sum(v(1:l,3))                                                  
      kct=0                                                             
      do 10901 i=1,l-1                                                  
      k=mt(i)                                                           
      sl=sl+v(k,2)                                                      
      sr=sr-v(k,2)                                                      
      wl=wl+v(k,3)                                                      
      wr=wr-v(k,3)                                                      
      if(wl.lt.fntn.or.wr.lt.fntn)goto 10901                            
      c=wr*wl*max(sr/wr,sl/wl)**pwr                                     
      if(c.le.cri)goto 10901                                            
      cri=c                                                             
      kct=i                                                             
      cri1=sl/wl                                                        
      cri2=sr/wr                                                        
      w1=wl                                                             
      w2=wr                                                             
      scri=max(cri1,cri2)                                               
10901 continue                                                          
      continue                                                          
      if(kct .ne. 0)goto 10921                                          
      cri=0.0                                                           
      return                                                            
10921 continue                                                          
      cat(1)=-kct                                                       
      do 10931 i=1,kct                                                  
      cat(i+1)=v(mt(i),1)                                               
10931 continue                                                          
      continue                                                          
      cri=scri                                                          
      kct=kct+1                                                         
      return                                                            
      end                                                               
      subroutine andarm(n,y,y2,z,w,dst,sw)                              
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      double precision y(n),y2(n),z(n),w(n)                             
      call set_kri(kri,2)                                               
      if(kri .ne. 1000)goto 10951                                       
      call rfcall(n, y, z, w, dst)                                      
      sw=sum(w)                                                         
      goto 10941                                                        
10951 if(kri .ne. 1)goto 10961                                          
      call andarm1(n,y,z,w,dst,sw)                                      
      goto 10941                                                        
10961 if(kri .ne. 2)goto 10971                                          
      call andarm2(n,y,z,w,dst,sw)                                      
      goto 10941                                                        
10971 if(kri .ne. 3)goto 10981                                          
      call andarm3(n,y,z,w,dst,sw)                                      
      goto 10941                                                        
10981 if(kri .ne. 4)goto 10991                                          
      call andarm4(n,y,z,w,dst,sw)                                      
      goto 10941                                                        
10991 if(kri .ne. 5)goto 11001                                          
      call andarm5(n,y,z,w,dst,sw)                                      
      goto 10941                                                        
11001 if(kri .ne. 6)goto 11011                                          
      call andarm6(n,y,y2,z,w,dst,sw)                                   
      goto 10941                                                        
11011 if(kri .ne. 7)goto 11021                                          
      call andarm7(n,y,z,w,dst,sw)                                      
      goto 10941                                                        
11021 if(kri .ne. 8)goto 11031                                          
      call andarm8(n,y,z,w,dst,sw)                                      
      goto 10941                                                        
11031 if(kri .ne. 9)goto 11041                                          
      call andarm7(n,y,z,w,dst,sw)                                      
      goto 10941                                                        
11041 if(kri .ne. 10)goto 11051                                         
      call andarm10(n,y,z,w,dst,sw)                                     
      goto 10941                                                        
11051 if(kri .ne. 11)goto 11061                                         
      call andarm11(dst,sw)                                             
      goto 10941                                                        
11061 if(kri .ne. 12)goto 11071                                         
      call andarm12(n,y,z,w,dst,sw)                                     
      goto 10941                                                        
11071 if(kri .ne. 13)goto 11081                                         
      call andarm12(n,y,z,w,dst,sw)                                     
      goto 10941                                                        
11081 if(kri .ne. 14)goto 11091                                         
      call andarm14(n,y,z,w,dst,sw)                                     
      goto 11101                                                        
11091 continue                                                          
      call andarm15(n,y,y2,z,w,dst,sw)                                  
11101 continue                                                          
10941 continue                                                          
      return                                                            
      end                                                               
      subroutine set_kri(irg,jrg)                                       
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      save kri                                                          
      if(jrg .ne. 1)goto 11121                                          
      kri=irg                                                           
      return                                                            
11121 continue                                                          
      irg=kri                                                           
      return                                                            
      end                                                               
      subroutine andarm11(dst,sw)                                       
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      dst=0.0                                                           
      sw=dst                                                            
      return                                                            
      end                                                               
      subroutine andarm2(n,y,z,w,dst,sw)                                
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      parameter(nmin=100)                                               
      double precision y(n),z(n),w(n)                                   
      integer my(n),mz(n)                                               
      call set_qqtrm(itrm,2)                                            
      if(n .ge. nmin)goto 11141                                         
      dst=0.0                                                           
      sw=sum(w)                                                         
      return                                                            
11141 continue                                                          
      if(n .ge. 2*itrm)goto 11161                                       
      dst=0.0                                                           
      sw=dst                                                            
      return                                                            
11161 continue                                                          
      do 11171 i=1,n                                                    
      my(i)=i                                                           
11171 continue                                                          
      continue                                                          
      call psort8(y,my,1,n)                                             
      do 11181 i=1,n                                                    
      mz(i)=i                                                           
11181 continue                                                          
      continue                                                          
      call psort8(z,mz,1,n)                                             
      dst=0.0                                                           
      sw1=dst                                                           
      do 11191 i=itrm+1,n-itrm                                          
      sw1=sw1+w(my(i))                                                  
      dst=dst+w(my(i))*abs(y(my(i))-z(mz(i)))                           
11191 continue                                                          
      continue                                                          
      dst=dst/sw1                                                       
      sw=sum(w)                                                         
      return                                                            
      end                                                               
      subroutine set_qqtrm(irg,jrg)                                     
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      save itrm                                                         
      if(jrg .ne. 1)goto 11211                                          
      itrm=irg                                                          
      return                                                            
11211 continue                                                          
      irg=itrm                                                          
      return                                                            
      end                                                               
      subroutine andarm1(n,y,z,w,dst,sw)                                
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      parameter(eps=1.0e-5,nmin=100)                                    
      double precision y(n),z(n),w(n),q(2*n),w2(2*n)                    
      integer m(2*n),iq(2*n)                                            
      if(n .ge. nmin)goto 11231                                         
      dst=0.0                                                           
      sw=sum(w)                                                         
      return                                                            
11231 continue                                                          
      n2=2*n                                                            
      do 11241 i=1,n                                                    
      q(i)=y(i)                                                         
      iq(i)=0                                                           
      q(i+n)=z(i)                                                       
      iq(i+n)=1                                                         
      w2(i)=w(i)                                                        
      w2(i+n)=w(i)                                                      
11241 continue                                                          
      continue                                                          
      do 11251 i=1,n2                                                   
      m(i)=i                                                            
11251 continue                                                          
      continue                                                          
      call psort8(q,m,1,n2)                                             
      sw=0.0                                                            
      tw=sw                                                             
      dst=tw                                                            
      sw2=2.0*sum(w)                                                    
      do 11261 i=1,n2                                                   
      k=m(i)                                                            
      if(iq(k) .ne. 0)goto 11281                                        
      sw=sw+w2(k)                                                       
      goto 11291                                                        
11281 continue                                                          
      tw=tw+w2(k)                                                       
11291 continue                                                          
      continue                                                          
      pw=(sw+tw)*(sw2-sw-tw)                                            
      pw=max(eps,pw)                                                    
      dst=dst+abs(sw-tw)/sqrt(pw)                                       
11261 continue                                                          
      continue                                                          
      dst=dst/n                                                         
      return                                                            
      end                                                               
      subroutine andarm3(n,y,z,w,dst,sw)                                
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      double precision y(n),z(n),w(n)                                   
      sw=sum(w)                                                         
      dst=dot_product(w,abs(y-z))/sw                                    
      return                                                            
      end                                                               
      subroutine andarm7(n,y,z,w,dst,sw)                                
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      parameter(nmin=20)                                                
      double precision y(n),z(n),w(n)                                   
      if(n .ge. nmin)goto 11311                                         
      dst=0.0                                                           
      sw=sum(w)                                                         
      return                                                            
11311 continue                                                          
      sw=sum(w)                                                         
      dst=abs(dot_product(w,y)/sw-dot_product(w,z)/sw)                  
      return                                                            
      end                                                               
      subroutine andarm12(n,y,z,w,dst,sw)                               
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      parameter(fmin=20)                                                
      double precision y(n),z(n),w(n)                                   
      if(n .ge. 2*int(fmin))goto 11331                                  
      dst=0.0                                                           
      sw=sum(w)                                                         
      return                                                            
11331 continue                                                          
      dst1=0.0                                                          
      dst2=dst1                                                         
      sw1=dst2                                                          
      sw2=sw1                                                           
      do 11341 i=1,n                                                    
      if(z(i) .ge. 0.0)goto 11361                                       
      sw1=sw1+w(i)                                                      
      dst1=dst1+w(i)*y(i)                                               
      goto 11371                                                        
11361 continue                                                          
      sw2=sw2+w(i)                                                      
      dst2=dst2+w(i)*y(i)                                               
11371 continue                                                          
      continue                                                          
11341 continue                                                          
      continue                                                          
      sw=sum(w)                                                         
      if((n*sw1/sw .ge. fmin) .and. (n*sw2/sw .ge. fmin))goto 11391     
      dst=0.0                                                           
      return                                                            
11391 continue                                                          
      dst=abs(dst2/sw2-dst1/sw1)                                        
      return                                                            
      end                                                               
      subroutine andarm14(n,y,z,w,dst,sw)                               
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      parameter(fmin=20,sml=-1.0e20)                                    
      double precision y(n),z(n),w(n)                                   
      if(n .ge. 2*int(fmin))goto 11411                                  
      dst=sml                                                           
      sw=sum(w)                                                         
      return                                                            
11411 continue                                                          
      dst1=0.0                                                          
      dst2=dst1                                                         
      sw1=dst2                                                          
      sw2=sw1                                                           
      do 11421 i=1,n                                                    
      if(z(i) .ge. 0.0)goto 11441                                       
      sw1=sw1+w(i)                                                      
      dst1=dst1+w(i)*y(i)                                               
      goto 11451                                                        
11441 continue                                                          
      sw2=sw2+w(i)                                                      
      dst2=dst2+w(i)*y(i)                                               
11451 continue                                                          
      continue                                                          
11421 continue                                                          
      continue                                                          
      sw=sum(w)                                                         
      if((n*sw1/sw .ge. fmin) .and. (n*sw2/sw .ge. fmin))goto 11471     
      dst=sml                                                           
      return                                                            
11471 continue                                                          
      dst=dst2/sw2-dst1/sw1                                             
      return                                                            
      end                                                               
      subroutine andarm8(n,y,z,w,dst,sw)                                
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      parameter(nmin=20,sml=-1.0e20)                                    
      double precision y(n),z(n),w(n)                                   
      if(n .ge. nmin)goto 11491                                         
      dst=sml                                                           
      sw=sum(w)                                                         
      return                                                            
11491 continue                                                          
      sw=sum(w)                                                         
      dst=dot_product(w,y)/sw-dot_product(w,z)/sw                       
      return                                                            
      end                                                               
      subroutine andarm4(n,y,z,w,dst,sw)                                
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      parameter(maxclass2=10000,nmin=100,idum=2)                        
      double precision y(n),z(n),w(n),out(maxclass2),dum(2,2)           
      double precision, dimension (:,:), allocatable :: costs           
      if(n .ge. nmin)goto 11511                                         
      dst=0.0                                                           
      sw=sum(w)                                                         
      return                                                            
11511 continue                                                          
      call classin(2,idum,dum,nclass,out)                               
      allocate(costs(1:nclass,1:nclass),stat=jerr);                     
      call reorg(2,nclass,out,costs)                                    
      dst=0.0                                                           
      do 11521 i=1,n                                                    
      ky=int(y(i)+0.1)                                                  
      kz=int(z(i)+0.1)                                                  
      dst=dst+w(i)*costs(ky,kz)                                         
11521 continue                                                          
      continue                                                          
      sw=sum(w)                                                         
      dst=dst/sw                                                        
      return                                                            
      end                                                               
      subroutine andarm5(n,y,z,w,dst,sw)                                
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      parameter(nmin=50)                                                
      double precision y(n),z(n),w(n)                                   
      data qntl /0.5/                                                   
      if(n .ge. nmin)goto 11541                                         
      dst=0.0                                                           
      sw=sum(w)                                                         
      return                                                            
11541 continue                                                          
      sw=sum(w)                                                         
      up=0.0                                                            
      do 11551 i=1,n                                                    
      if(y(i).le.z(i)) up=up+w(i)                                       
11551 continue                                                          
      continue                                                          
      dst=abs(up/sw-qntl)                                               
      return                                                            
      entry set_quant(arg)                                              
      qntl=arg                                                          
      return                                                            
      end                                                               
      subroutine andarm6(n,y,y2,z,w,dst,sw)                             
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      parameter(nit=100,xmiss=9.0e35,thr=1.0e-2,nmin=100)               
      double precision y(n),y2(n),z(n),w(n),yy(n,2)                     
      if(n .ge. nmin)goto 11571                                         
      dst=0.0                                                           
      sw=sum(w)                                                         
      return                                                            
11571 continue                                                          
      yy(:,1)=y                                                         
      yy(:,2)=y2                                                        
      call cendst(n,yy,z,w,nit,thr,xmiss,dst,sw)                        
      sw=sum(w)                                                         
      return                                                            
      end                                                               
      subroutine andarm15(n,y,y2,z,w,dst,sw)                            
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      parameter(nit=100,xmiss=9.0e35,thr=1.0e-2,nmin=100)               
      double precision y(n),y2(n),z(n),w(n),yy(n,2)                     
      if(n .ge. nmin)goto 11591                                         
      dst=0.0                                                           
      sw=sum(w)                                                         
      return                                                            
11591 continue                                                          
      yy(:,1)=y                                                         
      yy(:,2)=y2                                                        
      call cendst1(n,yy,z,w,nit,thr,xmiss,dst,sw)                       
      sw=sum(w)                                                         
      return                                                            
      end                                                               
      subroutine andarm10(n,y,z,w,dst,sw)                               
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      parameter(eps=1.0e-5,nmin=100)                                    
      double precision y(n),z(n),w(n)                                   
      integer m(n)                                                      
      sw=sum(w)                                                         
      if(n .ge. nmin)goto 11611                                         
      dst=0.0                                                           
      return                                                            
11611 continue                                                          
      sw1=0.0                                                           
      sw2=sw1                                                           
      do 11621 i=1,n                                                    
      m(i)=i                                                            
      if(z(i) .ge. 0.0)goto 11641                                       
      sw1=sw1+w(i)                                                      
      goto 11651                                                        
11641 continue                                                          
      sw2=sw2+w(i)                                                      
11651 continue                                                          
      continue                                                          
11621 continue                                                          
      continue                                                          
      call psort8(y,m,1,n)                                              
      s1=0.0                                                            
      s2=s1                                                             
      s=s2                                                              
      dst=s                                                             
      do 11661 i=1,n                                                    
      k=m(i)                                                            
      s=s+w(k)                                                          
      if(z(k) .ge. 0)goto 11681                                         
      s1=s1+w(k)/sw1                                                    
      goto 11691                                                        
11681 continue                                                          
      s2=s2+w(k)/sw2                                                    
11691 continue                                                          
      continue                                                          
      pw=s*(sw-s)                                                       
      pw=max(eps,pw)                                                    
      dst=dst+abs(s1-s2)/sqrt(pw)                                       
11661 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine stput (iseed)                                          
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      double precision x(*)                                             
      data i /987654321/                                                
      i=iseed                                                           
      return                                                            
      entry rget (x,n)                                                  
      do 1 j=1,n                                                        
c Naras fix: explicit conversion of nq to integer                       
c      i=mod(i*16807.0,2147483647.0)                                    
      i=int(mod(i*16807.0,2147483647.0))                                
      u=i                                                               
      u=u*.465661287d-9                                                 
c Naras fix: gcc fortran warns about label for statement in do loop     
c So put label on a separate continue statement                         
      x(j) = u                                                          
    1 continue                                                          
      return                                                            
      entry stget (irg)                                                 
      irg=i                                                             
      return                                                            
      end                                                               
      subroutine cendst(n,y,z,w,nit,thr,xmiss,dst,sw)                   
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      parameter(eps=0.1)                                                
      double precision y(n,2),z(n),w(n),b(2*n+1),q(3*n),cdf(3*n),r(n)   
      integer iq(3*n),mm(3*n),mz(n)                                     
      data nsamp /500/                                                  
      sw = sw + 0                                                       
      n2=2*n                                                            
      n3=3*n                                                            
      sw=sum(w)                                                         
      do 11701 i=1,n                                                    
      mz(i)=i                                                           
11701 continue                                                          
      continue                                                          
      call psort8(z,mz,1,n)                                             
      nq=int(0.25*n)                                                    
      teps=(z(mz(n-nq))-z(mz(nq)))*eps                                  
      do 11711 i=1,n                                                    
      if(y(i,2)-y(i,1).ge.teps)goto 11711                               
      y(i,1)=y(i,1)-teps                                                
      y(i,2)=y(i,2)+teps                                                
11711 continue                                                          
      continue                                                          
      do 11721 i=1,n                                                    
      b(i)=y(i,1)                                                       
      b(i+n)=y(i,2)                                                     
11721 continue                                                          
      continue                                                          
      m=0                                                               
      do 11731 i=1,n2                                                   
      if(b(i).le.-xmiss)goto 11731                                      
      if(b(i).ge.xmiss)goto 11731                                       
      m=m+1                                                             
      b(m)=b(i)                                                         
11731 continue                                                          
      continue                                                          
      call unique(m,b,nu)                                               
      if(nu .le. nsamp)goto 11751                                       
      call rget(r,nsamp)                                                
      do 11761 i=1,nsamp                                                
      r(i)=b(int(nu*r(i))+1)                                            
11761 continue                                                          
      continue                                                          
      nu=nsamp                                                          
      b(1:nu)=r(1:nu)                                                   
      call sort(b,nu)                                                   
11751 continue                                                          
      m=nu+1                                                            
      b(m)=xmiss                                                        
      call fintcdf1(n,y,m,b,w,nit,thr/m,cdf,jt,err)                     
      m=m-1                                                             
      do 11771 i=1,m                                                    
      q(i)=b(i)                                                         
      iq(i)=0                                                           
11771 continue                                                          
      continue                                                          
      do 11781 i=1,n                                                    
      mz(i)=i                                                           
11781 continue                                                          
      continue                                                          
      call psort8(z,mz,1,n)                                             
      mpn=m+n                                                           
      k=0                                                               
      do 11791 i=m+1,mpn                                                
      k=k+1                                                             
      q(i)=z(mz(k))                                                     
      iq(i)=1                                                           
11791 continue                                                          
      continue                                                          
      do 11801 i=1,mpn                                                  
      mm(i)=i                                                           
11801 continue                                                          
      continue                                                          
      call psort8(q,mm,1,mpn)                                           
      ycdf=0.0                                                          
      zcdf=ycdf                                                         
      dst=zcdf                                                          
      spw=dst                                                           
      ny=0                                                              
      nz=ny                                                             
      do 11811 i=1,mpn                                                  
      k=mm(i)                                                           
      if(iq(k) .ne. 0)goto 11831                                        
      ny=ny+1                                                           
      ycdf=cdf(ny)                                                      
      pw=float(i)*float(mpn-i)/float(mpn)**2                            
      pw=max(eps,pw)                                                    
      pw=1.0/sqrt(pw)                                                   
      spw=spw+pw                                                        
      dst=dst+pw*abs(ycdf-zcdf)                                         
      goto 11841                                                        
11831 continue                                                          
      nz=nz+1                                                           
      zcdf=zcdf+w(nz)/sw                                                
11841 continue                                                          
      continue                                                          
11811 continue                                                          
      continue                                                          
      dst=dst/spw                                                       
      return                                                            
      entry set_samp(irg)                                               
      nsamp=irg                                                         
      call set_samp1(nsamp)                                             
      return                                                            
      end                                                               
      subroutine cendst1(n,y,z,w,nit,thr,xmiss,dst,sw)                  
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      parameter(teps=0.01)                                              
      double precision y(n,2),z(n),w(n),b(2*n+1),cdf1(3*n),cdf2(3*n),r(n
     *)
      double precision y1(n,2),y2(n,2),w1(n),w2(n)                      
      data nsamp /500/                                                  
      n1=0.0                                                            
      n2=n1                                                             
      sw = sw + 0                                                       
      do 11851 i=1,n                                                    
      if(y(i,1).le.-xmiss)goto 11851                                    
      if(y(i,2).ge.xmiss)goto 11851                                     
      if(y(i,2)-y(i,1).ge.teps)goto 11851                               
      y(i,1)=y(i,1)-teps                                                
      y(i,2)=y(i,2)+teps                                                
11851 continue                                                          
      continue                                                          
      do 11861 i=1,n                                                    
      if(z(i) .ge. 0.0)goto 11881                                       
      n1=n1+1                                                           
      y1(n1,:)=y(i,:)                                                   
      w1(n1)=w(i)                                                       
      goto 11891                                                        
11881 continue                                                          
      n2=n2+1                                                           
      y2(n2,:)=y(i,:)                                                   
      w2(n2)=w(i)                                                       
11891 continue                                                          
      continue                                                          
11861 continue                                                          
      continue                                                          
      do 11901 i=1,n                                                    
      b(i)=y(i,1)                                                       
      b(i+n)=y(i,2)                                                     
11901 continue                                                          
      continue                                                          
      m=0                                                               
      do 11911 i=1,n2                                                   
      if(b(i).le.-xmiss)goto 11911                                      
      if(b(i).ge.xmiss)goto 11911                                       
      m=m+1                                                             
      b(m)=b(i)                                                         
11911 continue                                                          
      continue                                                          
      call unique(m,b,nu)                                               
      if(nu .le. nsamp)goto 11931                                       
      call rget(r,nsamp)                                                
      do 11941 i=1,nsamp                                                
      r(i)=b(int(nu*r(i))+1)                                            
11941 continue                                                          
      continue                                                          
      nu=nsamp                                                          
      b(1:nu)=r(1:nu)                                                   
      call sort(b,nu)                                                   
11931 continue                                                          
      m=nu+1                                                            
      b(m)=xmiss                                                        
      call getcdf1(n1,y1,w1,nit,thr,xmiss,nsamp,m,b,cdf1,sw1)           
      call getcdf1(n2,y2,w2,nit,thr,xmiss,nsamp,m,b,cdf2,sw2)           
      call diffcdf(m,cdf1,cdf2,dst)                                     
      return                                                            
      entry set_samp1(irg)                                              
      nsamp=irg                                                         
      return                                                            
      end                                                               
      subroutine getcdf1(n,y,w,nit,thr,xmiss,nsamp,m,b,cdf,sw)          
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      parameter(teps=0.01)                                              
      double precision y(n,2),w(n),b(2*n+1),cdf(3*n)                    
      xmiss = xmiss + 0                                                 
      nsamp = nsamp + 0                                                 
      n2=2*n                                                            
      sw=sum(w)                                                         
      call fintcdf1(n,y,m,b,w,nit,thr/m,cdf,jt,err)                     
      m=m-1                                                             
      return                                                            
      end                                                               
      subroutine diffcdf(m,cdf1,cdf2,dst)                               
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      double precision cdf1(m),cdf2(m)                                  
      f12=sqrt(float(m))                                                
      dst=0.0                                                           
      do 11951 i=1,m                                                    
      dst=dst+abs(cdf1(i)-cdf2(i))/sqrt(float(i)*float(m-i+1))          
11951 continue                                                          
      continue                                                          
      dst=f12*dst/m                                                     
      return                                                            
      end                                                               
      subroutine unique(n,y,nu)                                         
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      double precision y(n),yu(n)                                       
      integer m(n)                                                      
      do 11961 i=1,n                                                    
      m(i)=i                                                            
11961 continue                                                          
      continue                                                          
      call psort8(y,m,1,n)                                              
      nu=1                                                              
      yu(1)=y(m(1))                                                     
      do 11971 i=2,n                                                    
      if(y(m(i-1)).ge.y(m(i)))goto 11971                                
      nu=nu+1                                                           
      yu(nu)=y(m(i))                                                    
11971 continue                                                          
      continue                                                          
      y(1:nu)=yu(1:nu)                                                  
      return                                                            
      end                                                               
      subroutine fintcdf1(n,y,m,b,w1,nit,thr,cdf,jt,err)                
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      double precision y(n,2),b(m),w(n),w1(n),p(m),pij(n,m),ps(m),cdf(m)
      double precision djunk(1)                                         
      integer ijunk(1)                                                  
      integer, dimension (:), allocatable :: ic,jc,kc1,kc2,lc           
      double precision, dimension (:), allocatable :: smo               
      call set_vrb(ivrb,2)                                              
      w=w1/sum(w1)                                                      
      p=1.0/m                                                           
      nt=0                                                              
      do 11981 i=1,n                                                    
      do 11991 k=1,m                                                    
      if(y(i,1).ge.b(k))goto 11991                                      
      if(y(i,2).lt.b(k))goto 11991                                      
      nt=nt+1                                                           
11991 continue                                                          
      continue                                                          
11981 continue                                                          
      continue                                                          
      allocate(ic(1:(n+1)),stat=jerr)                                   
      allocate(jc(1:nt),stat=ierr)                                      
      jerr=jerr+ierr                                                    
      allocate(kc1(1:m),stat=ierr)                                      
      jerr=jerr+ierr                                                    
      allocate(kc2(1:m),stat=ierr)                                      
      jerr=jerr+ierr                                                    
      allocate(smo(1:m),stat=ierr)                                      
      jerr=jerr+ierr                                                    
      allocate(lc(1:nt),stat=ierr)                                      
      jerr=jerr+ierr                                                    
      if(jerr .eq. 0)goto 12011                                         
      err=8888.0                                                        
      return                                                            
12011 continue                                                          
      nt=0                                                              
      ic(1)=1                                                           
      do 12021 i=1,n                                                    
      do 12031 k=1,m                                                    
      if(y(i,1).ge.b(k))goto 12031                                      
      if(y(i,2).lt.b(k))goto 12031                                      
      nt=nt+1                                                           
      jc(nt)=k                                                          
12031 continue                                                          
      continue                                                          
      ic(i+1)=nt+1                                                      
12021 continue                                                          
      continue                                                          
      nt=0                                                              
      do 12041 j=1,m                                                    
      kc1(j)=nt+1                                                       
      do 12051 i=1,n                                                    
      if(y(i,1).ge.b(j))goto 12051                                      
      if(y(i,2).lt.b(j))goto 12051                                      
      nt=nt+1                                                           
      lc(nt)=i                                                          
12051 continue                                                          
      continue                                                          
      kc2(j)=nt                                                         
12041 continue                                                          
      continue                                                          
      if(ivrb .le. 0)goto 12071                                         
      call intpr('CDF iterations', -1, ijunk, 0)                        
12071 continue                                                          
      do 12081 it=1,nit                                                 
      jt=it                                                             
      ps=p                                                              
      do 12091 j=1,m                                                    
      pij(:,j)=0.0                                                      
      do 12101 ii=kc1(j),kc2(j)                                         
      i=lc(ii)                                                          
      s=sum(p(jc(ic(i):(ic(i+1)-1))))                                   
      if(s .gt. 0.0)goto 12121                                          
      err=-7777.0                                                       
      return                                                            
12121 continue                                                          
      pij(i,j)=w(i)*p(j)/s                                              
12101 continue                                                          
      continue                                                          
      p(j)=sum(pij(:,j))                                                
12091 continue                                                          
      continue                                                          
      if(m .le. 100)goto 12141                                          
      smo(1)=(2.0*p(1)+p(2))/3.0                                        
      smo(m)=(2.0*p(m)+p(m-1))/3.0                                      
      smo(2)=0.25*(p(1)+2.0*p(2)+p(3))                                  
      smo(m-1)=0.25*(p(m)+2.0*p(m-1)+p(m-2))                            
      do 12151 j=3,m-2                                                  
      smo(j)=(p(j-2)+2.0*p(j-1)+3.0*p(j)+2.0*p(j+1)+p(j+2))/9.0         
12151 continue                                                          
      continue                                                          
      p=smo                                                             
12141 continue                                                          
      err=sum(abs(p-ps))/m                                              
      if(kbad(err) .le. 0)goto 12171                                    
      err=7777.0                                                        
      return                                                            
12171 continue                                                          
      if(err.lt.thr)goto 12082                                          
      if(ivrb .le. 0)goto 12191                                         
      call intpr('.', -1, ijunk, 0)                                     
12191 continue                                                          
12081 continue                                                          
12082 continue                                                          
      cdf(1)=p(1)                                                       
      do 12201 j=2,m                                                    
      cdf(j)=cdf(j-1)+p(j)                                              
12201 continue                                                          
      continue                                                          
      if(ivrb .le. 0)goto 12221                                         
      djunk(1) = err                                                    
      call dblepr('Err = ', -1, djunk, 1)                               
12221 continue                                                          
      return                                                            
      end                                                               
      subroutine set_vrb(irg,jrg)                                       
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      save ivrb                                                         
      if(jrg .ne. 1)goto 12241                                          
      ivrb=irg                                                          
      return                                                            
12241 continue                                                          
      irg=ivrb                                                          
      return                                                            
      end                                                               
      subroutine cdfpoints1(m,x,n,y,w,cdf)                              
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      double precision x(m),y(n),w(n),cdf(m)                            
      i=1                                                               
      j=0                                                               
      sw=0.0                                                            
      continue                                                          
12251 continue                                                          
      j=j+1                                                             
      if(j.gt.m) go to 12260                                            
      continue                                                          
12271 if(y(i).gt.x(j))goto 12272                                        
      sw=sw+w(i)                                                        
      i=i+1                                                             
      if(i.gt.n)goto 12272                                              
      goto 12271                                                        
12272 continue                                                          
      if(i .le. n)goto 12291                                            
      do 12301 k=j,m                                                    
      cdf(k)=sw                                                         
12301 continue                                                          
      continue                                                          
      go to 12310                                                       
12291 continue                                                          
      cdf(j)=sw                                                         
      goto 12251                                                        
      continue                                                          
12310 continue                                                          
12260 continue                                                          
      cdf=cdf/sum(w)                                                    
      return                                                            
      end                                                               
      subroutine sort(x,n)                                              
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      double precision x(n),z(n)                                        
      integer m(n)                                                      
      do 12321 i=1,n                                                    
      m(i)=i                                                            
12321 continue                                                          
      continue                                                          
      z=x                                                               
      call psort8(z,m,1,n)                                              
      do 12331 i=1,n                                                    
      x(i)=z(m(i))                                                      
12331 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      function kbad(u)                                                  
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      kbad=0                                                            
      if(isnan(u).or.abs(u).ge.abs(huge(u))) kbad=1                     
      return                                                            
      end                                                               
      subroutine classin(ient,nclasssv,costssv,nout,out)                
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      double precision costssv(nclasssv,nclasssv),out(1)                
      double precision, dimension (:), allocatable :: costs             
      save costs,nclass                                                 
      nq=nclasssv*nclasssv                                              
      allocate(costs(1:nq),stat=jerr)                                   
      if(ient .ne. 1)goto 12351                                         
      nclass=nclasssv                                                   
      call reorg(1,nclass,costs,costssv)                                
      nout=1                                                            
      out=1.0                                                           
      goto 12361                                                        
12351 continue                                                          
      nout=nclass                                                       
      call reorg(2,nclass,costs,out)                                    
12361 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine reorg(ient,n,a,b)                                      
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      double precision a(n*n),b(n,n)                                    
      i=0                                                               
      if(ient .ne. 2)goto 12381                                         
      do 12391 k=1,n                                                    
      do 12401 j=1,n                                                    
      i=i+1                                                             
      b(j,k)=a(i)                                                       
12401 continue                                                          
      continue                                                          
12391 continue                                                          
      continue                                                          
      goto 12411                                                        
12381 continue                                                          
      do 12421 k=1,n                                                    
      do 12431 j=1,n                                                    
      i=i+1                                                             
      a(i)=b(j,k)                                                       
12431 continue                                                          
      continue                                                          
12421 continue                                                          
      continue                                                          
12411 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine crinode (itr,rtr,mxnodes,node,nodes,cri,wt)            
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      integer itr(6,*),nodes(mxnodes),m(mxnodes),ic(mxnodes)            
      double precision rtr(4,*),cri(mxnodes),wt(mxnodes),sc(mxnodes,2)  
      node=0                                                            
      k=itr(2,1)                                                        
      continue                                                          
12441 continue                                                          
      if(itr(4,k) .lt. 0)goto 12461                                     
      k=itr(2,k)                                                        
      goto 12441                                                        
12461 continue                                                          
      i1=itr(5,k)                                                       
      i2=itr(6,k)                                                       
      node=node+1                                                       
      if(node.gt.mxnodes) return                                        
      nodes(node)=k                                                     
      cri(node)=rtr(3,k)                                                
      wt(node)=rtr(4,k)                                                 
      continue                                                          
12471 if(k.eq.itr(2,iabs(itr(4,k))))goto 12472                          
      k=iabs(itr(4,k))                                                  
      if(k.eq.1)goto 12472                                              
      goto 12471                                                        
12472 continue                                                          
      if(k.eq.1)goto 12442                                              
      k=itr(3,iabs(itr(4,k)))                                           
      goto 12441                                                        
12442 continue                                                          
      do 12481 k=1,node                                                 
      m(k)=k                                                            
12481 continue                                                          
      continue                                                          
      call psort8(-cri,m,1,node)                                        
      do 12491 i=1,node                                                 
      ic(i)=nodes(m(i))                                                 
      sc(i,1)=cri(m(i))                                                 
      sc(i,2)=wt(m(i))                                                  
12491 continue                                                          
      continue                                                          
      do 12501 i=1,node                                                 
      nodes(i)=ic(i)                                                    
      cri(i)=sc(i,1)                                                    
      wt(i)=sc(i,2)                                                     
12501 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine prune1 (itr,rtr,nodes,thr,itro,rtro)                   
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      integer itr(6,nodes)                                              
      double precision rtr(4,nodes)                                     
      integer itro(6,nodes)                                             
      double precision rtro(4,nodes)                                    
      call prune(itr,rtr,nodes,thr)                                     
      itro=itr                                                          
      rtro=rtr                                                          
      return                                                            
      end                                                               
      subroutine prune (itr,rtr,nodes,thr)                              
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      integer itr(6,nodes)                                              
      double precision rtr(4,nodes)                                     
      continue                                                          
12511 continue                                                          
      nch=0                                                             
      do 12521 k=1,nodes                                                
      if(itr(4,k).le.0)goto 12521                                       
      nl=itr(2,k)                                                       
      nr=itr(3,k)                                                       
      if(itr(4,nl).ge.0)goto 12521                                      
      if(itr(4,nr).ge.0)goto 12521                                      
      if(max(rtr(3,nl),rtr(3,nr)).gt.rtr(3,k)+thr)goto 12521            
      itr(4,k)=-itr(4,k)                                                
      nch=nch+1                                                         
12521 continue                                                          
      continue                                                          
      if(nch.eq.0)goto 12512                                            
      goto 12511                                                        
12512 continue                                                          
      return                                                            
      end                                                               
      subroutine getnode (x,itre,rtre,cat,node)                         
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      integer itre(6,*)                                                 
      double precision x(*),rtre(4,*),cat(*)                            
      k=1                                                               
      continue                                                          
12531 if(itre(4,k).lt.0)goto 12532                                      
      if(itre(1,k) .le. 0)goto 12551                                    
      if(x(itre(1,k)) .ge. rtre(1,k))goto 12571                         
      k=itre(2,k)                                                       
      goto 12581                                                        
12571 continue                                                          
      k=itre(3,k)                                                       
12581 continue                                                          
      continue                                                          
      goto 12591                                                        
12551 continue                                                          
      j=-itre(1,k)                                                      
      kp=int(rtre(1,k)+0.1)                                             
      in=0                                                              
      np=int(abs(cat(kp))+0.1)                                          
      do 12601 i=1,np                                                   
      if(x(j).ne.cat(kp+i))goto 12601                                   
      in=1                                                              
      goto 12602                                                        
12601 continue                                                          
12602 continue                                                          
      if(cat(kp) .le. 0.0)goto 12621                                    
      if(in .ne. 0)goto 12641                                           
      k=itre(2,k)                                                       
      goto 12651                                                        
12641 continue                                                          
      k=itre(3,k)                                                       
12651 continue                                                          
      continue                                                          
      goto 12661                                                        
12621 continue                                                          
      if(in .eq. 0)goto 12681                                           
      k=itre(2,k)                                                       
      goto 12691                                                        
12681 continue                                                          
      k=itre(3,k)                                                       
12691 continue                                                          
      continue                                                          
12661 continue                                                          
      continue                                                          
12591 continue                                                          
      continue                                                          
      goto 12531                                                        
12532 continue                                                          
      node=k                                                            
      return                                                            
      end                                                               
      subroutine getnodes1 (no,ni,x,itre,rtre,cat,nodes)                
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      integer itre(6,*),nodes(no)                                       
      double precision x(no,ni),rtre(4,*),cat(*)                        
      do 12701 i=1,no                                                   
      call getnode (x(i,:),itre,rtre,cat,node)                          
      nodes(i)=node                                                     
12701 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine getlims(node,ni,itr,rtr,cat,nvar,jvar,vlims,jerr)      
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      integer itr(6,*),jvar(2,*)                                        
      double precision rtr(4,*),cat(*),vlims(*)                         
      ni = ni + 0                                                       
      jerr=0                                                            
      if(itr(4,node) .lt. 0)goto 12721                                  
      jerr=1                                                            
      return                                                            
12721 continue                                                          
      nvar=0                                                            
      k=node                                                            
      continue                                                          
12731 continue                                                          
      nvar=nvar+1                                                       
      kpp=abs(itr(4,k))                                                 
      if(itr(1,kpp) .le. 0)goto 12751                                   
      jvar(2,nvar)=0                                                    
      if(itr(2,kpp) .ne. k)goto 12771                                   
      jvar(1,nvar)=-itr(1,kpp)                                          
      goto 12781                                                        
12771 continue                                                          
      jvar(1,nvar)=itr(1,kpp)                                           
12781 continue                                                          
      continue                                                          
      vlims(nvar)=rtr(1,kpp)                                            
      goto 12791                                                        
12751 continue                                                          
      if(k .ne. itr(2,kpp))goto 12811                                   
      sgn=-1.0                                                          
      goto 12821                                                        
12811 continue                                                          
      sgn=1.0                                                           
12821 continue                                                          
      continue                                                          
      jvar(1,nvar)=-itr(1,kpp)                                          
      kp=int(rtr(1,kpp)+0.1)                                            
      jvar(2,nvar)=kp                                                   
      vlims(nvar)=sgn*abs(cat(kp))                                      
12791 continue                                                          
      continue                                                          
      k=kpp                                                             
      if(kpp.eq.1)goto 12732                                            
      goto 12731                                                        
12732 continue                                                          
      return                                                            
      end                                                               
      subroutine trans(ny,y,wy,nz,z,wz,nt,t)                            
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      double precision y(ny),wy(ny),z(nz),wz(nz),t(nt+2,2),u(max(ny,nz))
      double precision p(nt)                                            
      integer m(max(ny,nz))                                             
      do 12831 i=1,ny                                                   
      m(i)=i                                                            
      u(i)=y(i)                                                         
12831 continue                                                          
      continue                                                          
      call psort8(u,m,1,ny)                                             
      do 12841 i=1,ny                                                   
      y(i)=u(m(i))                                                      
12841 continue                                                          
      continue                                                          
      u=wy                                                              
      do 12851 i=1,ny                                                   
      wy(i)=u(m(i))                                                     
12851 continue                                                          
      continue                                                          
      do 12861 i=1,nz                                                   
      m(i)=i                                                            
      u(i)=z(i)                                                         
12861 continue                                                          
      continue                                                          
      call psort8(u,m,1,nz)                                             
      do 12871 i=1,nz                                                   
      z(i)=u(m(i))                                                      
12871 continue                                                          
      continue                                                          
      u=wz                                                              
      do 12881 i=1,nz                                                   
      wz(i)=u(m(i))                                                     
12881 continue                                                          
      continue                                                          
      do 12891 i=1,nt                                                   
      p(i)=(i-0.5)/float(nt)                                            
12891 continue                                                          
      continue                                                          
      call untie(ny,y,u)                                                
      call qntl(ny,u,wy,nt,p,t(:,1))                                    
      call untie(nz,z,u)                                                
      call qntl(nz,u,wz,nt,p,t(:,2))                                    
      return                                                            
      end                                                               
      subroutine qntl(n,y,w,nq,p,q)                                     
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      double precision y(n),w(n),p(nq),q(nq+2)                          
      sw=sum(w)                                                         
      k=1                                                               
      ff=w(1)                                                           
      q(1)=y(1)                                                         
      q(nq+2)=y(n)                                                      
      do 12901 i=2,n                                                    
      ff=ff+w(i)                                                        
      pp=ff/sw                                                          
      if(pp.lt.p(k))goto 12901                                          
      k=k+1                                                             
      q(k)=0.5*(y(i)+y(i-1))                                            
      if(k.ge.nq)goto 12902                                             
12901 continue                                                          
12902 continue                                                          
      q(nq+1)=0.5*(q(nq+2)+q(nq))                                       
      return                                                            
      end                                                               
      subroutine untie(n,y,u)                                           
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
      double precision y(n),u(n)                                        
      i=1                                                               
      k=0                                                               
      continue                                                          
12911 if(i.ge.n)goto 12912                                              
      if(y(i+1) .le. y(i))goto 12931                                    
      k=k+1                                                             
      u(k)=y(i)                                                         
      i=i+1                                                             
      goto 12911                                                        
12931 continue                                                          
      i1=i                                                              
      continue                                                          
12941 if(y(i+1).gt.y(i))goto 12942                                      
      i=i+1                                                             
      if(i.ge.n)goto 12942                                              
      goto 12941                                                        
12942 continue                                                          
      i2=i                                                              
      if(i1 .gt. 1)goto 12961                                           
      a=y(i1+1)                                                         
      b=y(i2+1)                                                         
      u(1)=y(1)                                                         
      k=1                                                               
      do 12971 j=i1+1,i2                                                
      k=k+1                                                             
      u(k)=a+(b-a)*(j-i1)/(i2-i1+1)                                     
12971 continue                                                          
      continue                                                          
      i=i2+1                                                            
      goto 12951                                                        
12961 if(i2 .lt. n)goto 12981                                           
      a=y(i1-1)                                                         
      b=(y(i2)-a)/(i2-i1+1)                                             
      do 12991 j=i1,i2                                                  
      k=k+1                                                             
      u(k)=a+b*(j-i1+1)                                                 
12991 continue                                                          
      continue                                                          
      goto 13001                                                        
12981 continue                                                          
      a=y(i1-1)                                                         
      b=y(i2)                                                           
      do 13011 j=i1,i2                                                  
      k=k+1                                                             
      u(k)=a+(b-a)*(j-i1+1)/(i2-i1+1)                                   
13011 continue                                                          
      continue                                                          
      i=i+1                                                             
13001 continue                                                          
12951 continue                                                          
      goto 12911                                                        
12912 continue                                                          
      if(k .ge. n)goto 13031                                            
      k=k+1                                                             
      u(k)=y(n)                                                         
13031 continue                                                          
      return                                                            
      end                                                               
      subroutine psort8 (v,a,ii,jj)                                     
      implicit double precision(a-h,o-z)                                
      implicit integer (i-n)
c     
c     puts into a the permutation vector which sorts v into             
c     increasing order. the array v is not modified.                    
c     only elements from ii to jj are considered.                       
c     arrays iu(k) and il(k) permit sorting up to 2**(k+1)-1 elements   
c                                                                       
c     this is a modification of cacm algorithm #347 by r. c. singleton, 
c     which is a modified hoare quicksort.                              
c                                                                       
      dimension a(jj),v(jj),iu(20),il(20)                               
      integer t,tt                                                      
      integer a                                                         
      double precision v                                                
      m=1                                                               
      i=ii                                                              
      j=jj                                                              
 10   if (i.ge.j) go to 80                                              
 20   k=i                                                               
      ij=(j+i)/2                                                        
      t=a(ij)                                                           
      vt=v(t)                                                           
      if (v(a(i)).le.vt) go to 30                                       
      a(ij)=a(i)                                                        
      a(i)=t                                                            
      t=a(ij)                                                           
      vt=v(t)                                                           
 30   l=j                                                               
      if (v(a(j)).ge.vt) go to 50                                       
      a(ij)=a(j)                                                        
      a(j)=t                                                            
      t=a(ij)                                                           
      vt=v(t)                                                           
      if (v(a(i)).le.vt) go to 50                                       
      a(ij)=a(i)                                                        
      a(i)=t                                                            
      t=a(ij)                                                           
      vt=v(t)                                                           
      go to 50                                                          
 40   a(l)=a(k)                                                         
      a(k)=tt                                                           
 50   l=l-1                                                             
      if (v(a(l)).gt.vt) go to 50                                       
      tt=a(l)                                                           
      vtt=v(tt)                                                         
 60   k=k+1                                                             
      if (v(a(k)).lt.vt) go to 60                                       
      if (k.le.l) go to 40                                              
      if (l-i.le.j-k) go to 70                                          
      il(m)=i                                                           
      iu(m)=l                                                           
      i=k                                                               
      m=m+1                                                             
      go to 90                                                          
 70   il(m)=k                                                           
      iu(m)=j                                                           
      j=l                                                               
      m=m+1                                                             
      go to 90                                                          
 80   m=m-1                                                             
      if (m.eq.0) return                                                
      i=il(m)                                                           
      j=iu(m)                                                           
 90   if (j-i.gt.10) go to 20                                           
      if (i.eq.ii) go to 10                                             
      i=i-1                                                             
 100  i=i+1                                                             
      if (i.eq.j) go to 80                                              
      t=a(i+1)                                                          
      vt=v(t)                                                           
      if (v(a(i)).le.vt) go to 100                                      
      k=i                                                               
 110  a(k+1)=a(k)                                                       
      k=k-1                                                             
      if (vt.lt.v(a(k))) go to 110                                      
      a(k+1)=t                                                          
      go to 100                                                         
      end                                                               
