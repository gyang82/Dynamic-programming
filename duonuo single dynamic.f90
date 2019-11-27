module input1
  implicit none
  type::tab01
    real::skH,krV
  end type tab01
  type::tab02
    real::xyH,ckQ
  end type tab02
  type::tab03
    real::syH,xxQ
  end type tab03
  !入库流量fq,12是它的列数，代表月份
  real,allocatable::fq(:,:)
  integer::N1,N2,N3,N4,NN4=12
contains
	 
integer function GetFileN( iFileUnit )
    Implicit None
    Integer , Intent( IN ) :: iFileUnit
    Integer :: ioS
    Character(Len=4) :: cDummy
    real :: cDummy1
	character ::filename
    GetFileN = 0	
    Rewind( iFileUnit )
    Do
      Read( iFileUnit ,*, ioStat = ioS ) cDummy1
	 ! filename=char(GetFileN)//'.txt'
	 ! open(unit=20,file='filename.txt',form='formatted')
	 ! write(20,*) cDummy1
      if ( ioS /= 0 ) Exit
      GetFileN = GetFileN + 1
    End Do
    Rewind( iFileUnit )
    Return
  End Function GetFileN
!这是一个线性插值函数，N为已有数据的长度，X，Y分别为已有的两组数据，已知T的值，对SS进行线型插值
  subroutine line(N,X,Y,T,SS)
	implicit none	
	integer N,i
	real X(N),Y(N),T,SS

	if(T==X(1)) then
		SS=Y(1);
	elseif(T==X(N))then
		SS=Y(N);
	else
		do i=1,N-1
			if((X(i)<=T.and.X(i+1)>=T).or.(X(i)>=T.and.X(i+1)<=T)) then
			SS=(Y(i+1)-Y(i))/(X(i+1)-X(i))*(T-X(i))+Y(i);
			end if
		end do
		if (T<X(1)) then
			SS=(Y(1+1)-Y(1))/(X(1+1)-X(1))*(T-X(1))+Y(1);
		elseif (T>X(N)) then
			SS=(Y(N)-Y(N-1))/(X(N)-X(N-1))*(T-X(N-1))+Y(N-1);
		end if
	end if
  end subroutine
  !效益函数，p为出力（kw），t为时间（s），bef为效益（kwh）
  real(kind=8) function bef(p,t,bzN)
	implicit none
	real p,t,bzN
	bef=p/1000
	if (p<bzN) then
		bef=-2*((bzN-p)/1000)**1.1 !100*((bzN-p)**1.8)发电量达到4.34出力有零     2*(money*(bzN-p)/1000)**1.1  
	end if
  end function bef
  subroutine maxx(n,a,maxa,b)
	implicit none 
	integer n,i,b
	real a(n),maxa
	maxa=maxval(a)
	do i=1,n
		if (a(i)==maxa) then
			b=i
			return
		end if
	end do
  end subroutine maxx
end
  
  
  
  
  program duonuo
  use input1
  implicit none
  !计算的时间段dt
  real,parameter::deadh=2320,normalh=2370,xunh=2366,xunh2=2368,xunh3=2370
  real::deadv,normalv,xunv,maxv
  real,allocatable:: h1(:,:),h2(:,:),qout(:,:),v1(:,:),v2(:,:),qin(:,:),clN(:,:),clNNN(:,:),qcap(:,:),mone(:,:),money(:,:),yy(:,:),maxh(:,:),hxia(:,:),hsha(:,:),qgo(:,:),qwai(:,:)
  !分别为出力系数，水头损失系数，时间（秒），单位换算系数，保证出力（kw），数据的列数
  real,parameter::k=8.69,ssh=0.02986,dt=24*30*3600,dw1=1e4,bzN=26000,zjN=100000,qmin(12)=1.8,basem=0.288
  integer,parameter::M4=12
  !y4用来统计保证率
  integer i,j,y4
  real dn1
  type(tab01),allocatable::Tab1(:)
  type(tab02),allocatable::Tab2(:)
  type(tab03),allocatable::Tab3(:)

!以下是由于动态规划法新添加的参数
  integer tt1,kk1,kk2,maxbef2n,bufzn
  real maxbef2
  integer,parameter::points1=100,points2=100
  real,allocatable:: maxf(:,:),cl1N(:,:),bef1(:,:),bef2(:,:),ddh1(:,:),bufh1(:),hh2(:,:)
  integer,allocatable::maxzn(:),maxfn(:,:)



  !call inputdata1()
  !读取水位库容关系曲线data01.txt
  open(5,file='data01.txt',status='old')
  N1=GetFileN(5)
  allocate(Tab1(N1))
  read(5,*) (Tab1(i)%skH,Tab1(i)%krV,i=1,N1)
  close(5)
  !读取尾水位和下泄流量关系曲线data02.txt
  open(5,file='data02.txt',status='old')
  N2=GetFileN(5)
  allocate(Tab2(N2))
  read(5,*) (Tab2(i)%xyH,Tab2(i)%ckQ,i=1,N2)
  close(5)
  !读取水位和下泄能力关系曲线data03.txt
  !在中长期调度中，由于时段长较长，受入、出库流量均化的影响，中长期发电调度无法考虑次洪水泄流控制，因此一般不需要考虑
  !泄流设施泄流能力的约束，这里输入的数据是假数据，目的是故意让泄流能力很大，达到不用考虑。
  open(5,file='data03.txt',status='old')
  N3=GetFileN(5)
  allocate(Tab3(N3))
  read(5,*) (Tab3(i)%syH,Tab3(i)%xxQ,i=1,N3)
  close(5)
  !读取96-05年流量资料data04.txt
  open(5,file='data04.txt',status='old')
  N4=GetFileN(5)
  allocate(fq(N4,NN4))
  do i=1,N4
    read(5,*) (fq(i,j),j=1,NN4)
  end do
  close(5)
  
  call line(N1,Tab1(:)%skH,Tab1(:)%krV,deadh,deadv);
  call line(N1,Tab1(:)%skH,Tab1(:)%krV,normalh,normalv);
  call line(N1,Tab1(:)%skH,Tab1(:)%krV,xunh,xunv);
  
  allocate(h1(N4,M4))
  allocate(h2(N4,M4))
  allocate(qout(N4,M4))
  allocate(v1(N4,M4))
  allocate(v2(N4,M4))
  allocate(qin(N4,M4))
  allocate(clN(N4,M4))
  allocate(clNNN(N4,M4))
  allocate(qcap(N4,M4))
  allocate(maxh(N4,M4))
  allocate(hxia(N4,M4))
  allocate(hsha(N4,M4))
  allocate(qgo(N4,M4))
  allocate(qwai(N4,M4))
   allocate(yy(N4,M4))
  allocate(money(N4,M4))
  allocate(mone(N4,M4))



  !由于动态规划增加的变量
  allocate(maxf(N4*M4+1,points1+1))
  allocate(maxfn(N4*M4+1,points1+1))
  allocate(cl1N(points1+1,points1+1))
  allocate(bef1(points1+1,points1+1))
  allocate(bef2(points1+1,points1+1))
  allocate(ddh1(N4,M4))
  allocate(bufh1(points1+1))
  allocate(hh2(N4*M4+1,points1+1))
  allocate(maxzn(N4*M4+1))

  !do i=1,N4
	  !h1(1,1)=0.5*(deadh+normalh)
	  h1(1,1)=2366
	  call line(N1,Tab1(:)%skH,Tab1(:)%krV,h1(1,1),v1(1,1))
  !end do
  !y4用来统计保证率
  y4=0
  !tt1用来统计状态的个数
  tt1=1
  do kk1=1,points1+1
	maxf(1,kk1)=0
	maxfn(1,kk1)=int(points1/2)
  end do
  do i=1,N4
	do j=1,M4
		qin(i,j)=fq(i,j)
		!这是一个辨别汛期的语句
			!这纴E尤丒说缂鄣挠跋丒
		if (mod(j,12)>=2.and.mod(j,12)<=6) then
			money(i,j)=basem*(1-0.24)
		else 
			if ((mod(j,12)>=1.and.mod(j,12)<=1).or.(mod(j,12)>=7.and.mod(j,12)<=7)) then
				money(i,j)=basem
			else
				money(i,j)=basem*(1+0.3)
			end if
		end if
		if (mod(j,12)>=1.and.mod(j,12)<=4) then
			if (mod(j,12)==4) then
				maxh(i,j)=xunh2
			!elseif (mod(j,12)==9) then
				!maxh(i,j)=xunh3
			else
				maxh(i,j)=xunh
			end if
		else
			maxh(i,j)=normalh
		end if
		ddh1(i,j)=(maxh(i,j)-deadh)/points1

		do kk1=1,points1+1
			bufh1(kk1)=deadh+(kk1-1)*ddh1(i,j)
		end do
		do kk1=1,points1+1
			hh2(tt1,kk1)=bufh1(kk1)
		end do
		do kk1=1,points1+1
			if (i==1.and.j==1) then
			!这个地方是动态规划的初始边界，也是时段初水库的起调水位
				!h1(i,j)=hh2(1,kk1)
				h1(i,j)=h1(i,j)
			else 
				h1(i,j)=hh2(tt1-1,kk1)
			end if
			call line(N1,Tab1(:)%skH,Tab1(:)%krV,h1(i,j),v1(i,j))
			do kk2=1,points1+1
				h2(i,j)=bufh1(kk2)		
				call line(N1,Tab1(:)%skH,Tab1(:)%krV,h2(i,j),v2(i,j))
				hsha(i,j)=0.5*(h1(i,j)+h2(i,j))
				call line(N3,Tab3(:)%syH,Tab3(:)%xxQ,hsha(i,j),qcap(i,j))
				qout(i,j)=(v1(i,j)-v2(i,j))*dw1/dt+qin(i,j)					
				call line(N2,Tab2(:)%ckQ,Tab2(:)%xyH,qout(i,j),hxia(i,j))
				cl1N(kk1,kk2)=k*qout(i,j)*(hsha(i,j)-hxia(i,j))*(1-ssh)
				if (cl1N(kk1,kk2)>zjN) then
					cl1N(kk1,kk2)=zjN
				elseif (cl1N(kk1,kk2)<0) then
					cl1N(kk1,kk2)=0
				end if
				
				bef1(kk1,kk2)=bef(cl1N(kk1,kk2),dt,bzN)
				
			!	if ((qout(i,j)<qmin(j)).or.(qout(i,j)>qcap(i,j))) then
				if ((qout(i,j)<qmin(j))) then
					!cl1N(kk1,kk2)=0
					bef1(kk1,kk2)=-1e8
				end if

				if (qout(i,j)>qcap(i,j).and.qcap(i,j)>=0) then
					bef1(kk1,kk2)=-1e8
				end if

				bef2(kk1,kk2)=maxf(tt1,kk1)+bef1(kk1,kk2)
			end do
		end do
		do kk2=1,points1+1
		!这是动态规划的初始条件，也是边界条件
			if (tt1==1) then
				maxf(tt1+1,kk2)=bef2(int(points1/2),kk2)
				maxfn(tt1+1,kk2)=int(points1/2)
			else
				call maxx(points1+1,bef2(:,kk2),maxbef2,maxbef2n)
				maxf(tt1+1,kk2)=maxbef2
				maxfn(tt1+1,kk2)=maxbef2n
			end if
		end do
		print *,tt1
		tt1=tt1+1
		
	end do 
  end do
  !这是动态规划的边界条件，也是调度的末水位，是很重要的控制条件
  call maxx(points1+1,maxf(tt1,:),maxbef2,maxbef2n)
  tt1=tt1-1;

  maxzn(tt1)=maxbef2n
  bufzn=maxbef2n
  do i=tt1-1,1,-1
	maxzn(i)=maxfn(i+2,bufzn)
	print *,maxfn(i+2,bufzn)
	bufzn=maxzn(i)
  end do
    print *, maxbef2
  tt1=1
  !h2(1,1)=deadh+(maxzn(tt1+1)-1)*ddh1(1,1)
  !call line(N1,Tab1(:)%skH,Tab1(:)%krV,h1(1,1),v1(1,1))
  h1(1,1)=2366
	  call line(N1,Tab1(:)%skH,Tab1(:)%krV,h1(1,1),v1(1,1))
  y4=0		
  do i=1,N4
	do j=1,M4
		if (tt1==1) then
			clNNN(i,j)=maxf(tt1+1,maxzn(tt1))-maxf(tt1,1)
		else
			clNNN(i,j)=maxf(tt1+1,maxzn(tt1))-maxf(tt1,maxzn(tt1-1))
		end if
		h2(i,j)=deadh+(maxzn(tt1)-1)*ddh1(i,j)
		call line(N1,Tab1(:)%skH,Tab1(:)%krV,h1(i,j),v1(i,j))
		call line(N1,Tab1(:)%skH,Tab1(:)%krV,h2(i,j),v2(i,j))
		qout(i,j)=(v1(i,j)-v2(i,j))*dw1/dt+qin(i,j)
		hsha(i,j)=0.5*(h1(i,j)+h2(i,j))
		call line(N2,Tab2(:)%ckQ,Tab2(:)%xyH,qout(i,j),hxia(i,j))
		clN(i,j)=k*qout(i,j)*(hsha(i,j)-hxia(i,j))*(1-ssh)
		if (clN(i,j)>zjN) then
			qgo(i,j)=zjN/(k*(hsha(i,j)-hxia(i,j))*(1-ssh))
			clN(i,j)=zjN
		elseif (clN(i,j)<0) then
			qgo(i,j)=0
			clN(i,j)=0
		else 
			qgo(i,j)=qout(i,j)
		end if
		qwai(i,j)=qout(i,j)-qgo(i,j)
		!万元
		mone(i,j)=clN(i,j)*money(i,j)*dt/3600/1e4
		if (clN(i,j)-bzN>-0.5) then
			y4=y4+1
			yy(i,j)=1
		else
			yy(i,j)=0
		end if
		if (j<M4) then
			v1(i,j+1)=v2(i,j)
			h1(i,j+1)=h2(i,j)
		else
			if (i<N4) then
				v1(i+1,1)=v2(i,j)
				h1(i+1,1)=h2(i,j)
			end if
		end if
		!这是回代过程的点睛之笔，杨光曾经因为这颗“钉子”彻夜难眠，后人一定要谨记，否则后果可能更严重！！！
		tt1=tt1+1
	end do
  end do
	print*,'以下是利用动态规划法进行调度的调度结果：'
	print*,"如果保证出力为26000kW，则保证率为：",real(y4*1./(1.*N4*M4));
	print*,'年发电量为：',sum(clN(:,:))/N4*dt/3600/1e8,'亿千瓦时';
	print*,'7-12月发电量为：',sum(clN(:,3:8))/N4*dt/3600/1e8,'亿千瓦时';
	print*,'年发电量为：',sum(clN(:,:))/1000,'亿千瓦时';
	print*,'年发电效益为',sum(mone(:,:))/N4,'万元';
	print*,'合计为：',sum(clN(:,:))*dt/3600/1e8,'亿千瓦时';
	print*,'年弃水量为：',sum(qwai(:,:))/N4*dt/1e8,'亿立方米';
	print*,'7-12月弃水量为：',sum(qwai(:,3:8))/N4*dt/1e8,'亿立方米';

	open(2001,file='clNd01.txt');
	open(2002,file='H1d01.txt');
	open(2003,file='H2d01.txt');
	open(2004,file='qoutd01.txt');
	open(2005,file='qgod01.txt');
	open(2006,file='v201.txt');
	open(2007,file='v101.txt');
	open(2008,file='yyd01.txt');

	do i=1,N4
		write(2001,'(1x,f10.3\)')(clN(i,j),j=1,M4);
		write(2001,'(/\)')
	end do
	do i=1,N4
		write(2002,'(1x,f10.3\)')(h1(i,j),j=1,M4)
		write(2002,'(/\)')
	end do
	do i=1,N4
		write(2003,'(1x,f10.3\)')(h2(i,j),j=1,M4)
		write(2003,'(/\)')
	end do

	do i=1,N4
		write(2004,'(1x,f10.3\)')(qout(i,j),j=1,M4)
		write(2004,'(/\)')
	end do;
	do i=1,N4
		write(2005,'(1x,f10.3\)')(qgo(i,j),j=1,M4)
		write(2005,'(/\)')
	end do;
	do i=1,N4
		write(2006,'(1x,f10.3\)')(v2(i,j),j=1,M4)
		write(2006,'(/\)')
	end do;
	do i=1,N4
		write(2007,'(1x,f10.3\)')(v1(i,j),j=1,M4)
		write(2007,'(/\)')
	end do;
	do i=1,N4
		write(2008,'(1x,f10.3\)')(yy(i,j),j=1,M4)
		write(2008,'(/\)')
	end do;
	!写入文件
	!这里值得注意的是怎样让它输入的数据不换行
	!因为fortran在输入数据时会默认在一定长度后换行,这是为了方便阅读，但有时候我们并不需要它
	!	pause 
  !dt=90
	
   
end

