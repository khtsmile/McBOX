	  o1  w   k820309    h          14.0        ŕŠë\                                                                                                           
       physics.f90 PHYSICS                   @                               
                                                          
                            @                              
                                                           
                         @                               
                                                           
       KEFF K_COL K_TL          @                                         
       FISSION_BANK THREAD_BANK BANK_IDX                   @                               'X             &      #ID 	   #N_COORD 
   #CELL_INSTANCE    #COORD    #LAST_N_COORD    #LAST_CELL    #E    #LAST_E    #G    #LAST_G    #WGT     #MU !   #ALIVE "   #LAST_XYZ_CURRENT #   #LAST_XYZ $   #LAST_UVW %   #LAST_WGT &   #ABSORB_WGT '   #FISSION (   #EVENT )   #EVENT_NUCLIDE *   #EVENT_MT +   #DELAYED_GROUP ,   #N_BANK -   #WGT_BANK .   #N_DELAYED_BANK /   #SURFACE 0   #CELL_BORN 1   #MATERIAL 2   #LAST_MATERIAL 3   #SQRTKT 4   #LAST_SQRTKT 5   #N_COLLISION 6   #WRITE_TRACK 7   #YES_SAB 8   #CLEAR 9   #INITIALIZE <   #SET ?                 $                             	                                 $                              
                                $                                                              $                                   
              P             #LOCALCOORD    p          p 
           p 
                                        @  @                              'P              
      #CELL    #UNIVERSE    #LATTICE    #LATTICE_X    #LATTICE_Y    #LATTICE_Z    #XYZ    #UVW    #DIST    #RESET                                                                                                                                                0                                                                                                                                               0                                                                                                                                               0                                                                                                                                               0                                                                                                                                               0                                                                                                                                               0                                                                            
  p          p            p                                                                                 0                 
  p          p            p                                                                          H       	   
   1         Ŕ                                             
     #RESET_COORD    #         H     @                                               #THIS              
                                     P               #LOCALCOORD                  $                                   0                          $                                   
       4                  p          p 
           p 
                                       $                                  `         
                 $                                  h         
                 $                                   p      	                    $                                   t      
                    $                                   x         
                 $                             !              
                 $                              "                               $                             #                            
  p          p            p                                        $                             $            ¨                
  p          p            p                                        $                             %            Ŕ                
  p          p            p                                        $                             &     Ř         
                 $                             '     ŕ         
                 $                              (     č                          $                              )     ě                          $                              *     đ                          $                              +     ô                          $                              ,     ř                          $                              -     ü                          $                             .               
                 $                              /                              p          p            p                                        $                              0     (                          $                              1     ,                          $                              2     0                          $                              3     4                          $                             4     8         
                 $                             5     @          
                 $                              6     H      !                   $                              7     L      "                                                                                                      $                              8     P      #                                                                                         1         Ŕ    $                            9             $     #CLEAR :   #         @     @                            :                    #THIS ;                                             ;     X              #PARTICLE    1         Ŕ    $                            <             %     #INITIALIZE =   #         @     @                            =                    #THIS >                                             >     X              #PARTICLE    1         Ŕ    $                            ?             &     #SET_PARTICLE @   #         @     @                            @                    #THIS A   #SOURCE B                                             A     X              #PARTICLE                                               B     H               #BANK C                  @  @                          D     'Ř                    #BASE E   #LEN F   #OFFSET G   #FLAGS H   #RANK I   #RESERVED1 J   #DIMINFO K                                              E                                                              F                                                             G                                                             H                                                             I                                                              J     (                                                         K            0                    #FOR_DESC_TRIPLET L   p          p            p                                         @  @                         L     '                    #EXTENT M   #MULT N   #LOWERBOUND O                                              M                                                              N                                                             O                                    @                          P     '                                       @  @                          C     'H                    #WGT Q   #XYZ R   #UVW S   #E T   #G U                                              Q                
                                              R                             
  p          p            p                                                                     S                              
  p          p            p                                                                     T     8          
                                               U     @                               @               @           V     '(                   #MAT_ID W   #SIG_TR X   #SIG_ABS Y   #SIG_CAP Z   #SIG_FIS [   #NU \   #CHI ]   #SIG_SCAT ^                                              W                                                                    X                             
            &                                                                                    Y            `                 
            &                                                                                    Z            ¨                 
            &                                                                                    [            đ                 
            &                                                                                    \            8                
            &                                                                                   ]                            
            &                                                                                   ^            Č                
            &                   &                                                                                       _     
                  @                                `     
                                                   a     
                                                  b            H                        &                                           #BANK C              @                                c     č      H              p          p č          p č                        #BANK C              @                                d                                                       e            (                       &                                           #MG_XS V   %         @                               f                     
       n                         V              Crang                    (         `                                g                                   
    #RAND_VEC%SQRT h   #RAND_VEC%COS i   #RAND_VEC%SIN j   p          p            p                                                                      h     SQRT                                             i     COS                                             j     SIN                                             k     
                
                 {ŽGáz?        1E-2#         @                                   l                   #COLLISION_MG%SIZE m   #COLLISION_MG%INT n   #COLLISION_MG%SUM o   #P p                                              m     SIZE                                            n     INT                                            o     SUM           
D                                 p     X              #PARTICLE    #         @                                   q                   #COLLISION_MG_DT%SIZE r   #COLLISION_MG_DT%INT s   #COLLISION_MG_DT%SUM t   #P u   #MACRO_MAJOR v                                              r     SIZE                                            s     INT                                            t     SUM           
D                                 u     X              #PARTICLE              
                                 v     
                   fn#fn    ź   @   J   OMP_LIB    ü   @   J   CONSTANTS     <  @   J   PARTICLE_HEADER    |  @   J   XS_HEADER    ź  @   J   RANDOMS    ü  P   J  VARIABLES    L  b   J  BANK_HEADER )   Ž  ]      PARTICLE+PARTICLE_HEADER ,     H   a   PARTICLE%ID+PARTICLE_HEADER 1   S  H   a   PARTICLE%N_COORD+PARTICLE_HEADER 7     H   a   PARTICLE%CELL_INSTANCE+PARTICLE_HEADER /   ă  Ź   a   PARTICLE%COORD+PARTICLE_HEADER +     É      LOCALCOORD+PARTICLE_HEADER 0   X  Ľ   a   LOCALCOORD%CELL+PARTICLE_HEADER 4   ý  Ľ   a   LOCALCOORD%UNIVERSE+PARTICLE_HEADER 3   ˘  Ľ   a   LOCALCOORD%LATTICE+PARTICLE_HEADER 5   G	  Ľ   a   LOCALCOORD%LATTICE_X+PARTICLE_HEADER 5   ě	  Ľ   a   LOCALCOORD%LATTICE_Y+PARTICLE_HEADER 5   
  Ľ   a   LOCALCOORD%LATTICE_Z+PARTICLE_HEADER /   6     a   LOCALCOORD%XYZ+PARTICLE_HEADER /   Ň     a   LOCALCOORD%UVW+PARTICLE_HEADER 0   n  H   a   LOCALCOORD%DIST+PARTICLE_HEADER 1   ś  Y   a   LOCALCOORD%RESET+PARTICLE_HEADER ,     R      RESET_COORD+PARTICLE_HEADER 1   a  X   a   RESET_COORD%THIS+PARTICLE_HEADER 6   š  H   a   PARTICLE%LAST_N_COORD+PARTICLE_HEADER 3        a   PARTICLE%LAST_CELL+PARTICLE_HEADER +     H   a   PARTICLE%E+PARTICLE_HEADER 0   ĺ  H   a   PARTICLE%LAST_E+PARTICLE_HEADER +   -  H   a   PARTICLE%G+PARTICLE_HEADER 0   u  H   a   PARTICLE%LAST_G+PARTICLE_HEADER -   ˝  H   a   PARTICLE%WGT+PARTICLE_HEADER ,     H   a   PARTICLE%MU+PARTICLE_HEADER /   M  H   a   PARTICLE%ALIVE+PARTICLE_HEADER :        a   PARTICLE%LAST_XYZ_CURRENT+PARTICLE_HEADER 2   1     a   PARTICLE%LAST_XYZ+PARTICLE_HEADER 2   Í     a   PARTICLE%LAST_UVW+PARTICLE_HEADER 2   i  H   a   PARTICLE%LAST_WGT+PARTICLE_HEADER 4   ą  H   a   PARTICLE%ABSORB_WGT+PARTICLE_HEADER 1   ů  H   a   PARTICLE%FISSION+PARTICLE_HEADER /   A  H   a   PARTICLE%EVENT+PARTICLE_HEADER 7     H   a   PARTICLE%EVENT_NUCLIDE+PARTICLE_HEADER 2   Ń  H   a   PARTICLE%EVENT_MT+PARTICLE_HEADER 7     H   a   PARTICLE%DELAYED_GROUP+PARTICLE_HEADER 0   a  H   a   PARTICLE%N_BANK+PARTICLE_HEADER 2   Š  H   a   PARTICLE%WGT_BANK+PARTICLE_HEADER 8   ń     a   PARTICLE%N_DELAYED_BANK+PARTICLE_HEADER 1     H   a   PARTICLE%SURFACE+PARTICLE_HEADER 3   Ő  H   a   PARTICLE%CELL_BORN+PARTICLE_HEADER 2     H   a   PARTICLE%MATERIAL+PARTICLE_HEADER 7   e  H   a   PARTICLE%LAST_MATERIAL+PARTICLE_HEADER 0   ­  H   a   PARTICLE%SQRTKT+PARTICLE_HEADER 5   ő  H   a   PARTICLE%LAST_SQRTKT+PARTICLE_HEADER 5   =  H   a   PARTICLE%N_COLLISION+PARTICLE_HEADER 5     ¤   a   PARTICLE%WRITE_TRACK+PARTICLE_HEADER 1   )  ¤   a   PARTICLE%YES_SAB+PARTICLE_HEADER /   Í  S   a   PARTICLE%CLEAR+PARTICLE_HEADER &      R      CLEAR+PARTICLE_HEADER +   r  V   a   CLEAR%THIS+PARTICLE_HEADER 4   Č  X   a   PARTICLE%INITIALIZE+PARTICLE_HEADER +      R      INITIALIZE+PARTICLE_HEADER 0   r  V   a   INITIALIZE%THIS+PARTICLE_HEADER -   Č  Z   a   PARTICLE%SET+PARTICLE_HEADER -   "  ^      SET_PARTICLE+PARTICLE_HEADER 2     V   a   SET_PARTICLE%THIS+PARTICLE_HEADER 4   Ö  R   a   SET_PARTICLE%SOURCE+PARTICLE_HEADER 3   (         FOR_ARRAY_DESCRIPTOR+ISO_C_BINDING 8   Č  H   a   FOR_ARRAY_DESCRIPTOR%BASE+ISO_C_BINDING 7     H   a   FOR_ARRAY_DESCRIPTOR%LEN+ISO_C_BINDING :   X  H   a   FOR_ARRAY_DESCRIPTOR%OFFSET+ISO_C_BINDING 9      H   a   FOR_ARRAY_DESCRIPTOR%FLAGS+ISO_C_BINDING 8   č  H   a   FOR_ARRAY_DESCRIPTOR%RANK+ISO_C_BINDING =   0  H   a   FOR_ARRAY_DESCRIPTOR%RESERVED1+ISO_C_BINDING ;   x  ˛   a   FOR_ARRAY_DESCRIPTOR%DIMINFO+ISO_C_BINDING /   *  v      FOR_DESC_TRIPLET+ISO_C_BINDING 6      H   a   FOR_DESC_TRIPLET%EXTENT+ISO_C_BINDING 4   č  H   a   FOR_DESC_TRIPLET%MULT+ISO_C_BINDING :   0   H   a   FOR_DESC_TRIPLET%LOWERBOUND+ISO_C_BINDING '   x   P       #UNLPOLY+ISO_C_BINDING !   Č   y      BANK+BANK_HEADER %   A!  H   a   BANK%WGT+BANK_HEADER %   !     a   BANK%XYZ+BANK_HEADER %   %"     a   BANK%UVW+BANK_HEADER #   Á"  H   a   BANK%E+BANK_HEADER #   	#  H   a   BANK%G+BANK_HEADER     Q#  Ž       MG_XS+XS_HEADER '   ˙#  P   a   MG_XS%MAT_ID+XS_HEADER '   O$     a   MG_XS%SIG_TR+XS_HEADER (   ă$     a   MG_XS%SIG_ABS+XS_HEADER (   w%     a   MG_XS%SIG_CAP+XS_HEADER (   &     a   MG_XS%SIG_FIS+XS_HEADER #   &     a   MG_XS%NU+XS_HEADER $   3'     a   MG_XS%CHI+XS_HEADER )   Ç'  Ź   a   MG_XS%SIG_SCAT+XS_HEADER    s(  @       KEFF+VARIABLES     ł(  @       K_COL+VARIABLES    ó(  @       K_TL+VARIABLES )   3)         FISSION_BANK+BANK_HEADER (   É)         THREAD_BANK+BANK_HEADER %   g*  @       BANK_IDX+BANK_HEADER     §*         XS_MG+XS_HEADER    >+         RANG+RANDOMS !   ×+  Ű       RAND_VEC+RANDOMS &   ˛,  =      RAND_VEC%SQRT+RANDOMS %   ď,  <      RAND_VEC%COS+RANDOMS %   +-  <      RAND_VEC%SIN+RANDOMS "   g-  t       WGT_MIN+CONSTANTS    Ű-         COLLISION_MG "   m.  =      COLLISION_MG%SIZE !   Ş.  <      COLLISION_MG%INT !   ć.  <      COLLISION_MG%SUM    "/  V   a   COLLISION_MG%P     x/  Ź       COLLISION_MG_DT %   $0  =      COLLISION_MG_DT%SIZE $   a0  <      COLLISION_MG_DT%INT $   0  <      COLLISION_MG_DT%SUM "   Ů0  V   a   COLLISION_MG_DT%P ,   /1  @   a   COLLISION_MG_DT%MACRO_MAJOR 