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
  Integer Function GetFileN( iFileUnit )
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
	  open(unit=20,file='filename.txt',form='formatted')
	  write(20,*) cDummy1
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
end