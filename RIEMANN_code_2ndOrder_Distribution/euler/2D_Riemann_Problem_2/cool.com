! ----------------------------------------------------------------------

! This is "cool.com"

      INTEGER ncool

      PARAMETER ( ncool = 500)

      REAL tabtemp, tempmin,                                            &
     &     truetemp ( ncool), coolfn ( ncool), dcoolfn ( ncool),        &
     &     heatfn ( ncool), dheatfn ( ncool)

      COMMON /COOLCOM/ tabtemp, tempmin,                                &
     &       truetemp, coolfn, dcoolfn,                                 &
     &       heatfn, dheatfn

! ----------------------------------------------------------------------


