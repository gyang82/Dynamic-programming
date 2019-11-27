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