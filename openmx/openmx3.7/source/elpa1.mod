  ©O  Å   k820309    4          12.1        yÛQ                                                                                                           
       elpa1.f90 ELPA1              TRIDIAG_REAL TRANS_EV_REAL MULT_AT_B_REAL TRIDIAG_COMPLEX TRANS_EV_COMPLEX MULT_AH_B_COMPLEX SOLVE_TRIDI CHOLESKY_REAL INVERT_TRM_REAL CHOLESKY_COMPLEX INVERT_TRM_COMPLEX LOCAL_INDEX LEAST_COMMON_MULTIPLE HH_TRANSFORM_REAL HH_TRANSFORM_COMPLEX TIME_EVP_FWD TIME_EVP_SOLVE TIME_EVP_BACK ELPA_PRINT_TIMES #         @                                                 	   #TRIDIAG_REAL%MIN    #TRIDIAG_REAL%DOT_PRODUCT    #TRIDIAG_REAL%UBOUND    #TRIDIAG_REAL%MAX    #TRIDIAG_REAL%MOD    #NA    #A    #LDA 	   #NBLK 
   #MPI_COMM_ROWS    #MPI_COMM_COLS    #D 
   #E    #TAU                  @                                 MIN               @                                 DOT_PRODUCT               @                                 UBOUND               @                                 MAX               @                                 MOD           D @@                                                 B  D @@                                                 
       p        5  p        r 	   p          5  p        r 	     1     5  p        r 	     1                             D @@                              	                      D @@                              
                      D @@                                                    D @@                                                   D @@                             
                    
     p          5  p        r        5  p        r                               D @@                                                 
     p          5  p        r        5  p        r                               D @@                                                 
     p          5  p        r        5  p        r                      #         @                                                 
   #TRANS_EV_REAL%MIN    #TRANS_EV_REAL%UBOUND    #TRANS_EV_REAL%MAX    #TRANS_EV_REAL%MOD    #NA    #NQC    #A    #LDA    #TAU    #Q    #LDQ    #NBLK    #MPI_COMM_ROWS    #MPI_COMM_COLS                  @                                 MIN               @                                 UBOUND               @                                 MAX               @                                 MOD            @@                                                    D @@                                                 B    @                                                 
       p        5  p        r    p          5  p        r      1     5  p        r      1                               @                                                     @                                                 
     p          5  p        r        5  p        r                             B  D @@                                                 
       p        5  p        r    p          5  p        r      1     5  p        r      1                             D @@                                                    D @@                                                    D @@                                                    D @@                                          #         @                                                 
   #MULT_AT_B_REAL%MIN     #MULT_AT_B_REAL%UBOUND !   #MULT_AT_B_REAL%MOD "   #UPLO_A #   #UPLO_C $   #NA %   #NCB &   #A '   #LDA (   #B )   #LDB *   #NBLK +   #MPI_COMM_ROWS ,   #MPI_COMM_COLS -   #C .   #LDC /                 @                                  MIN               @                            !     UBOUND               @                            "     MOD             @                              #                                        @                              $                                      D @@                              %                      D @@                              &                   B    @                             '                    
       p        5  p        r (   p          5  p        r (     1     5  p        r (     1                               @                              (                   B  D @@                             )                    
       p        5  p        r *   p          5  p        r *     1     5  p        r *     1                             D @@                              *                      D @@                              +                      D @@                              ,                      D @@                              -                   B  D  @                             .                    
       p        5  p 
       r /   p          5  p 
       r /     1     5  p 
       r /     1                               @                              /            #         @                                  0               	   #TRIDIAG_COMPLEX%CONJG 1   #TRIDIAG_COMPLEX%MIN 2   #TRIDIAG_COMPLEX%DOT_PRODUCT 3   #TRIDIAG_COMPLEX%UBOUND 4   #TRIDIAG_COMPLEX%MAX 5   #TRIDIAG_COMPLEX%MOD 6   #NA 7   #A 8   #LDA 9   #NBLK :   #MPI_COMM_ROWS ;   #MPI_COMM_COLS <   #D =   #E >   #TAU ?                 @                            1     CONJG               @                            2     MIN               @                            3     DOT_PRODUCT               @                            4     UBOUND               @                            5     MAX               @                            6     MOD           D @@                              7                   B  D @@                             8                    
 %      p        5  p        r 9   p          5  p        r 9     1     5  p        r 9     1                             D @@                              9                      D @@                              :                      D @@                              ;                      D @@                              <                     D @@                             =                    
 '    p          5  p        r 7       5  p        r 7                              D @@                             >                    
 (    p          5  p        r 7       5  p        r 7                              D @@                             ?                    
   &    p          5  p        r 7       5  p        r 7                     #         @                                  @               
     #TRANS_EV_COMPLEX%CONJG A   #TRANS_EV_COMPLEX%MIN B   #TRANS_EV_COMPLEX%UBOUND C   #TRANS_EV_COMPLEX%MAX D   #TRANS_EV_COMPLEX%MOD E   #NA F   #NQC G   #A H   #LDA I   #TAU J   #Q K   #LDQ L   #NBLK M   #MPI_COMM_ROWS N   #MPI_COMM_COLS O                 @                            A     CONJG               @                            B     MIN               @                            C     UBOUND               @                            D     MAX               @                            E     MOD            @@                              F                      D @@                              G                   B    @                             H                    
   4      p        5  p        r I   p          5  p        r I     1     5  p        r I     1                               @                              I                       @                             J                    
   6    p          5  p        r F       5  p        r F                            B  D @@                             K                    
   5      p        5  p        r L   p          5  p        r L     1     5  p        r L     1                             D @@                              L                      D @@                              M                      D @@                              N                      D @@                              O            #         @                                  P               
     #MULT_AH_B_COMPLEX%MIN Q   #MULT_AH_B_COMPLEX%UBOUND R   #MULT_AH_B_COMPLEX%MOD S   #UPLO_A T   #UPLO_C U   #NA V   #NCB W   #A X   #LDA Y   #B Z   #LDB [   #NBLK \   #MPI_COMM_ROWS ]   #MPI_COMM_COLS ^   #C _   #LDC `                 @                            Q     MIN               @                            R     UBOUND               @                            S     MOD             @                              T                                        @                              U                                      D @@                              V                      D @@                              W                   B    @                             X                    
   >      p        5  p        r Y   p          5  p        r Y     1     5  p        r Y     1                               @                              Y                   B  D @@                             Z                    
   ?      p        5  p        r [   p          5  p        r [     1     5  p        r [     1                             D @@                              [                      D @@                              \                      D @@                              ]                      D @@                              ^                   B  D  @                             _                    
           @      p        5  p 
                 r `   p          5  p 
                 r `     1     5  p 
                 r `     1                               @                              `            #         @                                  a               	   #SOLVE_TRIDI%ABS b   #SOLVE_TRIDI%MIN c   #NA d   #NEV e   #D f   #E g   #Q h   #LDQ i   #NBLK j   #MPI_COMM_ROWS k   #MPI_COMM_COLS l                 @                            b     ABS               @                            c     MIN           D @@                              d                       @@                              e                     D @@                             f                    
           G    p          5  p        r d       5  p        r d                              D @@                             g                    
           H    p          5  p        r d       5  p        r d                            B  D @@                             h                    
           I      p        5  p        r i   p          5  p        r i     1     5  p        r i     1                             D @@                              i                      D @@                              j                      D @@                              k                      D @@                              l            #         @                                  m                  #CHOLESKY_REAL%MIN n   #CHOLESKY_REAL%UBOUND o   #CHOLESKY_REAL%MAX p   #CHOLESKY_REAL%MOD q   #NA r   #A s   #LDA t   #NBLK u   #MPI_COMM_ROWS v   #MPI_COMM_COLS w                 @                            n     MIN               @                            o     UBOUND               @                            p     MAX               @                            q     MOD           D @@                              r                   B  D @@                             s                    
                 p        5  p        r t   p          5  p        r t     1     5  p        r t     1                             D @@                              t                      D @@                              u                      D @@                              v                      D @@                              w            #         @                                  x                  #INVERT_TRM_REAL%UBOUND y   #INVERT_TRM_REAL%MOD z   #NA {   #A |   #LDA }   #NBLK ~   #MPI_COMM_ROWS    #MPI_COMM_COLS                  @                            y     UBOUND               @                            z     MOD           D @@                              {                   B  D @@                             |                    
                         p        5  p        r 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     }   p          5  p        r 
}     1     5  p        r 
}     1                             D @@                              
}                      D @@                              ~                      D @@                                                    D @@                                          #         @                                                    #CHOLESKY_COMPLEX%CONJG    #CHOLESKY_COMPLEX%MIN    #CHOLESKY_COMPLEX%UBOUND    #CHOLESKY_COMPLEX%MAX    #CHOLESKY_COMPLEX%MOD    #NA    #A    #LDA    #NBLK    #MPI_COMM_ROWS    #MPI_COMM_COLS                  @                                 CONJG               @                                 MIN               @                                 UBOUND               @                                 MAX               @                                 MOD           D @@                                                 B  D @@                                                 
         p        5  p        r    p          5  p        r      1     5  p        r      1                             D @@                                                    D @@                                                    D @@                                                    D @@                                          #         @                                                    #INVERT_TRM_COMPLEX%UBOUND    #INVERT_TRM_COMPLEX%MOD    #NA    #A    #LDA    #NBLK    #MPI_COMM_ROWS    #MPI_COMM_COLS                  @                                 UBOUND               @                                 MOD           D @@                                                 B  D @@                                                 
         p        5  p        r    p          5  p        r      1     5  p        r      1                             D @@                                                    D @@                                                    D @@                                                    D @@                                          %         @                                    
                       #LOCAL_INDEX%MOD    #IDX    #MY_PROC    #NUM_PROCS    #NBLK    #IFLAG                  @                                 MOD             @                                                      @                                                     @@                                                     @@                                                      @                                          %         @                                    
                       #LEAST_COMMON_MULTIPLE%MOD    #A    #B                   @                                 MOD           
     @                                                   
    @@                                          #         @                                 ¡                  #HH_TRANSFORM_REAL%SQRT ¢   #HH_TRANSFORM_REAL%SIGN £   #ALPHA ¤   #XNORM_SQ ¥   #XF ¦   #TAU §                 @                            ¢     SQRT               @                            £     SIGN           
  D @@                             ¤     
                   
     @                             ¥     
                  D  @                             ¦     
                   D  @                             §     
         #         @                                 ¨                  #HH_TRANSFORM_COMPLEX%DCMPLX ©   #HH_TRANSFORM_COMPLEX%DIMAG ª   #HH_TRANSFORM_COMPLEX%DBLE «   #HH_TRANSFORM_COMPLEX%SQRT ¬   #HH_TRANSFORM_COMPLEX%SIGN ­   #ALPHA ®   #XNORM_SQ ¯   #XF °   #TAU ±                 @                            ©     DCMPLX               @                            ª     DIMAG               @                            «     DBLE               @                            ¬     SQRT               @                            ­     SIGN           
  D @@                             ®     
                   
     @                             ¯     
                  D  @                             °     
                   D  @                             ±     
         
                                              ²     
         
                                              ³     
         
                                              ´     
         
                                             µ                           @                           ¶                          #MPI_BOTTOM ·                @                           ·                                   @                           ¸                          #MPI_IN_PLACE ¹                @                           ¹                                   @                           º                          #MPI_ARGV_NULL »   -             @                           »                                 p          p            p                                                 @                           ¼                          #MPI_ARGVS_NULL ½                @                           ½             
                        @                           ¾                          #MPI_ERRCODES_IGNORE ¿                @                           ¿                                 p          p            p                                                 @                           À                          #MPI_STATUS_IGNORE Á                @                           Á                                 p          p            p                                                 @                           Â                          #MPI_STATUSES_IGNORE Ã                @                           Ã             
                      fn#fn    ¸   ?  b   uapp(ELPA1    ÷         TRIDIAG_REAL !     <      TRIDIAG_REAL%MIN )   S  D      TRIDIAG_REAL%DOT_PRODUCT $     ?      TRIDIAG_REAL%UBOUND !   Ö  <      TRIDIAG_REAL%MAX !     <      TRIDIAG_REAL%MOD     N  @   a   TRIDIAG_REAL%NA      ô   a   TRIDIAG_REAL%A !     @   a   TRIDIAG_REAL%LDA "   Â  @   a   TRIDIAG_REAL%NBLK +     @   a   TRIDIAG_REAL%MPI_COMM_ROWS +   B  @   a   TRIDIAG_REAL%MPI_COMM_COLS      ´   a   TRIDIAG_REAL%D    6  ´   a   TRIDIAG_REAL%E !   ê  ´   a   TRIDIAG_REAL%TAU            TRANS_EV_REAL "   ¯	  <      TRANS_EV_REAL%MIN %   ë	  ?      TRANS_EV_REAL%UBOUND "   *
    <      TRANS_EV_REAL%MAX "   f
    <      TRANS_EV_REAL%MOD !   ¢
    @   a   TRANS_EV_REAL%NA "   â
    @   a   TRANS_EV_REAL%NQC     "  ô   a   TRANS_EV_REAL%A "     @   a   TRANS_EV_REAL%LDA "   V  ´   a   TRANS_EV_REAL%TAU     
    ô   a   TRANS_EV_REAL%Q "   þ
    @   a   TRANS_EV_REAL%LDQ #   >  @   a   TRANS_EV_REAL%NBLK ,   ~  @   a   TRANS_EV_REAL%MPI_COMM_ROWS ,   ¾  @   a   TRANS_EV_REAL%MPI_COMM_COLS    þ        MULT_AT_B_REAL #     <      MULT_AT_B_REAL%MIN &   V  ?      MULT_AT_B_REAL%UBOUND #     <      MULT_AT_B_REAL%MOD &   Ñ  P   a   MULT_AT_B_REAL%UPLO_A &   !  P   a   MULT_AT_B_REAL%UPLO_C "   q  @   a   MULT_AT_B_REAL%NA #   ±  @   a   MULT_AT_B_REAL%NCB !   ñ  ô   a   MULT_AT_B_REAL%A #   å  @   a   MULT_AT_B_REAL%LDA !   %  ô   a   MULT_AT_B_REAL%B #     @   a   MULT_AT_B_REAL%LDB $   Y  @   a   MULT_AT_B_REAL%NBLK -     @   a   MULT_AT_B_REAL%MPI_COMM_ROWS -   Ù  @   a   MULT_AT_B_REAL%MPI_COMM_COLS !     ô   a   MULT_AT_B_REAL%C #   
    @   a   MULT_AT_B_REAL%LDC     M  J      TRIDIAG_COMPLEX &     >      TRIDIAG_COMPLEX%CONJG $   Õ  <      TRIDIAG_COMPLEX%MIN ,     D      TRIDIAG_COMPLEX%DOT_PRODUCT '   U  ?      TRIDIAG_COMPLEX%UBOUND $     <      TRIDIAG_COMPLEX%MAX $   Ð  <      TRIDIAG_COMPLEX%MOD #     @   a   TRIDIAG_COMPLEX%NA "   L  ô   a   TRIDIAG_COMPLEX%A $   @  @   a   TRIDIAG_COMPLEX%LDA %     @   a   TRIDIAG_COMPLEX%NBLK .   À  @   a   TRIDIAG_COMPLEX%MPI_COMM_ROWS .      @   a   TRIDIAG_COMPLEX%MPI_COMM_COLS "   @  ´   a   TRIDIAG_COMPLEX%D "   ô  ´   a   TRIDIAG_COMPLEX%E $   ¨  ´   a   TRIDIAG_COMPLEX%TAU !   \  9      TRANS_EV_COMPLEX '     >      TRANS_EV_COMPLEX%CONJG %   Ó  <      TRANS_EV_COMPLEX%MIN (     ?      TRANS_EV_COMPLEX%UBOUND %   N  <      TRANS_EV_COMPLEX%MAX %     <      TRANS_EV_COMPLEX%MOD $   Æ  @   a   TRANS_EV_COMPLEX%NA %      @   a   TRANS_EV_COMPLEX%NQC #   F   ô   a   TRANS_EV_COMPLEX%A %   :!  @   a   TRANS_EV_COMPLEX%LDA %   z!  ´   a   TRANS_EV_COMPLEX%TAU #   ."  ô   a   TRANS_EV_COMPLEX%Q %   "#  @   a   TRANS_EV_COMPLEX%LDQ &   b#  @   a   TRANS_EV_COMPLEX%NBLK /   ¢#  @   a   TRANS_EV_COMPLEX%MPI_COMM_ROWS /   â#  @   a   TRANS_EV_COMPLEX%MPI_COMM_COLS "   "$  %      MULT_AH_B_COMPLEX &   G%  <      MULT_AH_B_COMPLEX%MIN )   %  ?      MULT_AH_B_COMPLEX%UBOUND &   Â%  <      MULT_AH_B_COMPLEX%MOD )   þ%  P   a   MULT_AH_B_COMPLEX%UPLO_A )   N&  P   a   MULT_AH_B_COMPLEX%UPLO_C %   &  @   a   MULT_AH_B_COMPLEX%NA &   Þ&  @   a   MULT_AH_B_COMPLEX%NCB $   '  ô   a   MULT_AH_B_COMPLEX%A &   (  @   a   MULT_AH_B_COMPLEX%LDA $   R(  ô   a   MULT_AH_B_COMPLEX%B &   F)  @   a   MULT_AH_B_COMPLEX%LDB '   )  @   a   MULT_AH_B_COMPLEX%NBLK 0   Æ)  @   a   MULT_AH_B_COMPLEX%MPI_COMM_ROWS 0   *  @   a   MULT_AH_B_COMPLEX%MPI_COMM_COLS $   F*  ô   a   MULT_AH_B_COMPLEX%C &   :+  @   a   MULT_AH_B_COMPLEX%LDC    z+  Ñ       SOLVE_TRIDI     K,  <      SOLVE_TRIDI%ABS     ,  <      SOLVE_TRIDI%MIN    Ã,  @   a   SOLVE_TRIDI%NA     -  @   a   SOLVE_TRIDI%NEV    C-  ´   a   SOLVE_TRIDI%D    ÷-  ´   a   SOLVE_TRIDI%E    «.  ô   a   SOLVE_TRIDI%Q     /  @   a   SOLVE_TRIDI%LDQ !   ß/  @   a   SOLVE_TRIDI%NBLK *   0  @   a   SOLVE_TRIDI%MPI_COMM_ROWS *   _0  @   a   SOLVE_TRIDI%MPI_COMM_COLS    0  ï       CHOLESKY_REAL "   1  <      CHOLESKY_REAL%MIN %   Ê1  ?      CHOLESKY_REAL%UBOUND "   	2  <      CHOLESKY_REAL%MAX "   E2  <      CHOLESKY_REAL%MOD !   2  @   a   CHOLESKY_REAL%NA     Á2  ô   a   CHOLESKY_REAL%A "   µ3  @   a   CHOLESKY_REAL%LDA #   õ3  @   a   CHOLESKY_REAL%NBLK ,   54  @   a   CHOLESKY_REAL%MPI_COMM_ROWS ,   u4  @   a   CHOLESKY_REAL%MPI_COMM_COLS     µ4  Å       INVERT_TRM_REAL '   z5  ?      INVERT_TRM_REAL%UBOUND $   ¹5  <      INVERT_TRM_REAL%MOD #   õ5  @   a   INVERT_TRM_REAL%NA "   56  ô   a   INVERT_TRM_REAL%A $   )7  @   a   INVERT_TRM_REAL%LDA %   i7  @   a   INVERT_TRM_REAL%NBLK .   ©7  @   a   INVERT_TRM_REAL%MPI_COMM_ROWS .   é7  @   a   INVERT_TRM_REAL%MPI_COMM_COLS !   )8        CHOLESKY_COMPLEX '   @9  >      CHOLESKY_COMPLEX%CONJG %   ~9  <      CHOLESKY_COMPLEX%MIN (   º9  ?      CHOLESKY_COMPLEX%UBOUND %   ù9  <      CHOLESKY_COMPLEX%MAX %   5:  <      CHOLESKY_COMPLEX%MOD $   q:  @   a   CHOLESKY_COMPLEX%NA #   ±:  ô   a   CHOLESKY_COMPLEX%A %   ¥;  @   a   CHOLESKY_COMPLEX%LDA &   å;  @   a   CHOLESKY_COMPLEX%NBLK /   %<  @   a   CHOLESKY_COMPLEX%MPI_COMM_ROWS /   e<  @   a   CHOLESKY_COMPLEX%MPI_COMM_COLS #   ¥<  Ë       INVERT_TRM_COMPLEX *   p=  ?      INVERT_TRM_COMPLEX%UBOUND '   ¯=  <      INVERT_TRM_COMPLEX%MOD &   ë=  @   a   INVERT_TRM_COMPLEX%NA %   +>  ô   a   INVERT_TRM_COMPLEX%A '   ?  @   a   INVERT_TRM_COMPLEX%LDA (   _?  @   a   INVERT_TRM_COMPLEX%NBLK 1   ?  @   a   INVERT_TRM_COMPLEX%MPI_COMM_ROWS 1   ß?  @   a   INVERT_TRM_COMPLEX%MPI_COMM_COLS    @         LOCAL_INDEX     ¾@  <      LOCAL_INDEX%MOD     ú@  @   a   LOCAL_INDEX%IDX $   :A  @   a   LOCAL_INDEX%MY_PROC &   zA  @   a   LOCAL_INDEX%NUM_PROCS !   ºA  @   a   LOCAL_INDEX%NBLK "   úA  @   a   LOCAL_INDEX%IFLAG &   :B  }       LEAST_COMMON_MULTIPLE *   ·B  <      LEAST_COMMON_MULTIPLE%MOD (   óB  @   a   LEAST_COMMON_MULTIPLE%A (   3C  @   a   LEAST_COMMON_MULTIPLE%B "   sC  ª       HH_TRANSFORM_REAL '   D  =      HH_TRANSFORM_REAL%SQRT '   ZD  =      HH_TRANSFORM_REAL%SIGN (   D  @   a   HH_TRANSFORM_REAL%ALPHA +   ×D  @   a   HH_TRANSFORM_REAL%XNORM_SQ %   E  @   a   HH_TRANSFORM_REAL%XF &   WE  @   a   HH_TRANSFORM_REAL%TAU %   E        HH_TRANSFORM_COMPLEX ,   §F  ?      HH_TRANSFORM_COMPLEX%DCMPLX +   æF  >      HH_TRANSFORM_COMPLEX%DIMAG *   $G  =      HH_TRANSFORM_COMPLEX%DBLE *   aG  =      HH_TRANSFORM_COMPLEX%SQRT *   G  =      HH_TRANSFORM_COMPLEX%SIGN +   ÛG  @   a   HH_TRANSFORM_COMPLEX%ALPHA .   H  @   a   HH_TRANSFORM_COMPLEX%XNORM_SQ (   [H  @   a   HH_TRANSFORM_COMPLEX%XF )   H  @   a   HH_TRANSFORM_COMPLEX%TAU    ÛH  @       TIME_EVP_FWD    I  @       TIME_EVP_SOLVE    [I  @       TIME_EVP_BACK !   I  @       ELPA_PRINT_TIMES )   ÛI  `      ELPA1!MPI_FORTRAN_BOTTOM    ;J  H      MPI_BOTTOM +   J  b      ELPA1!MPI_FORTRAN_IN_PLACE    åJ  H      MPI_IN_PLACE ,   -K  c      ELPA1!MPI_FORTRAN_ARGV_NULL    K  ¤      MPI_ARGV_NULL -   4L  d      ELPA1!MPI_FORTRAN_ARGVS_NULL    L  H      MPI_ARGVS_NULL 2   àL  i      ELPA1!MPI_FORTRAN_ERRCODES_IGNORE $   IM  ¤      MPI_ERRCODES_IGNORE 0   íM  g      ELPA1!MPI_FORTRAN_STATUS_IGNORE "   TN  ¤      MPI_STATUS_IGNORE 2   øN  i      ELPA1!MPI_FORTRAN_STATUSES_IGNORE $   aO  H      MPI_STATUSES_IGNORE 