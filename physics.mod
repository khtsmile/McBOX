	  à1  x   k820309    h          14.0        Ví]                                                                                                           
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
      #CELL    #UNIVERSE    #LATTICE    #LATTICE_X    #LATTICE_Y    #LATTICE_Z    #XYZ    #UVW    #DIST    #RESET                                                                                                                                                0                                                                                                                                               0                                                                                                                                               0                                                                                                                                               0                                                                                                                                               0                                                                                                                                               0                                                                            
  p          p            p                                                                                 0                 
  p          p            p                                                                          H       	   
   1         À                                             
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
  p          p            p                                        $                             %            À                
  p          p            p                                        $                             &     Ø         
                 $                             '     à         
                 $                              (     è                          $                              )     ì                          $                              *     ð                          $                              +     ô                          $                              ,     ø                          $                              -     ü                          $                             .               
                 $                              /                              p          p            p                                        $                              0     (                          $                              1     ,                          $                              2     0                          $                              3     4                          $                             4     8         
                 $                             5     @          
                 $                              6     H      !                   $                              7     L      "                                                                                                      $                              8     P      #                                                                                         1         À    $                            9             $     #CLEAR :   #         @     @                            :                    #THIS ;                                             ;     X              #PARTICLE    1         À    $                            <             %     #INITIALIZE =   #         @     @                            =                    #THIS >                                             >     X              #PARTICLE    1         À    $                            ?             &     #SET_PARTICLE @   #         @     @                            @                    #THIS A   #SOURCE B                                             A     X              #PARTICLE                                               B     H               #BANK C                  @  @                          D     'Ø                    #BASE E   #LEN F   #OFFSET G   #FLAGS H   #RANK I   #RESERVED1 J   #DIMINFO K                                              E                                                              F                                                             G                                                             H                                                             I                                                              J     (                                                         K            0                    #FOR_DESC_TRIPLET L   p          p            p                                         @  @                         L     '                    #EXTENT M   #MULT N   #LOWERBOUND O                                              M                                                              N                                                             O                                    @                          P     '                                       @  @                          C     'H                    #WGT Q   #XYZ R   #UVW S   #E T   #G U                                              Q                
                                              R                             
  p          p            p                                                                     S                              
  p          p            p                                                                     T     8          
                                               U     @                               @               @           V     '(                   #MAT_ID W   #SIG_TR X   #SIG_ABS Y   #SIG_CAP Z   #SIG_FIS [   #NU \   #CHI ]   #SIG_SCAT ^                                              W                                                                    X                             
            &                                                                                    Y            `                 
            &                                                                                    Z            ¨                 
            &                                                                                    [            ð                 
            &                                                                                    \            8                
            &                                                                                   ]                            
            &                                                                                   ^            È                
            &                   &                                                                                       _     
                  @                                `     
                                                   a     
                                                  b            H                        &                                           #BANK C              @                                c     è      H              p          p è          p è                        #BANK C              @                                d                                                       e            (                       &                                           #MG_XS V   %         @                               f                     
       n                         V              Crang                    (         `                                g                                   
    #RAND_VEC%SQRT h   #RAND_VEC%COS i   #RAND_VEC%SIN j   p          p            p                                                                      h     SQRT                                             i     COS                                             j     SIN                                             k     
                
                 {®Gáz?        1E-2                                             l                                                       0#         @                                   m                   #COLLISION_MG%SIZE n   #COLLISION_MG%INT o   #COLLISION_MG%SUM p   #P q                                              n     SIZE                                            o     INT                                            p     SUM           
D                                 q     X              #PARTICLE    #         @                                   r                   #COLLISION_MG_DT%SIZE s   #COLLISION_MG_DT%INT t   #COLLISION_MG_DT%SUM u   #P v   #MACRO_MAJOR w                                              s     SIZE                                            t     INT                                            u     SUM           
D                                 v     X              #PARTICLE              
                                 w     
                   fn#fn    ¼   @   J   OMP_LIB    ü   @   J   CONSTANTS     <  @   J   PARTICLE_HEADER    |  @   J   XS_HEADER    ¼  @   J   RANDOMS    ü  P   J  VARIABLES    L  b   J  BANK_HEADER )   ®  ]      PARTICLE+PARTICLE_HEADER ,     H   a   PARTICLE%ID+PARTICLE_HEADER 1   S  H   a   PARTICLE%N_COORD+PARTICLE_HEADER 7     H   a   PARTICLE%CELL_INSTANCE+PARTICLE_HEADER /   ã  ¬   a   PARTICLE%COORD+PARTICLE_HEADER +     É      LOCALCOORD+PARTICLE_HEADER 0   X  ¥   a   LOCALCOORD%CELL+PARTICLE_HEADER 4   ý  ¥   a   LOCALCOORD%UNIVERSE+PARTICLE_HEADER 3   ¢  ¥   a   LOCALCOORD%LATTICE+PARTICLE_HEADER 5   G	  ¥   a   LOCALCOORD%LATTICE_X+PARTICLE_HEADER 5   ì	  ¥   a   LOCALCOORD%LATTICE_Y+PARTICLE_HEADER 5   
  ¥   a   LOCALCOORD%LATTICE_Z+PARTICLE_HEADER /   6     a   LOCALCOORD%XYZ+PARTICLE_HEADER /   Ò     a   LOCALCOORD%UVW+PARTICLE_HEADER 0   n  H   a   LOCALCOORD%DIST+PARTICLE_HEADER 1   ¶  Y   a   LOCALCOORD%RESET+PARTICLE_HEADER ,     R      RESET_COORD+PARTICLE_HEADER 1   a  X   a   RESET_COORD%THIS+PARTICLE_HEADER 6   ¹  H   a   PARTICLE%LAST_N_COORD+PARTICLE_HEADER 3        a   PARTICLE%LAST_CELL+PARTICLE_HEADER +     H   a   PARTICLE%E+PARTICLE_HEADER 0   å  H   a   PARTICLE%LAST_E+PARTICLE_HEADER +   -  H   a   PARTICLE%G+PARTICLE_HEADER 0   u  H   a   PARTICLE%LAST_G+PARTICLE_HEADER -   ½  H   a   PARTICLE%WGT+PARTICLE_HEADER ,     H   a   PARTICLE%MU+PARTICLE_HEADER /   M  H   a   PARTICLE%ALIVE+PARTICLE_HEADER :        a   PARTICLE%LAST_XYZ_CURRENT+PARTICLE_HEADER 2   1     a   PARTICLE%LAST_XYZ+PARTICLE_HEADER 2   Í     a   PARTICLE%LAST_UVW+PARTICLE_HEADER 2   i  H   a   PARTICLE%LAST_WGT+PARTICLE_HEADER 4   ±  H   a   PARTICLE%ABSORB_WGT+PARTICLE_HEADER 1   ù  H   a   PARTICLE%FISSION+PARTICLE_HEADER /   A  H   a   PARTICLE%EVENT+PARTICLE_HEADER 7     H   a   PARTICLE%EVENT_NUCLIDE+PARTICLE_HEADER 2   Ñ  H   a   PARTICLE%EVENT_MT+PARTICLE_HEADER 7     H   a   PARTICLE%DELAYED_GROUP+PARTICLE_HEADER 0   a  H   a   PARTICLE%N_BANK+PARTICLE_HEADER 2   ©  H   a   PARTICLE%WGT_BANK+PARTICLE_HEADER 8   ñ     a   PARTICLE%N_DELAYED_BANK+PARTICLE_HEADER 1     H   a   PARTICLE%SURFACE+PARTICLE_HEADER 3   Õ  H   a   PARTICLE%CELL_BORN+PARTICLE_HEADER 2     H   a   PARTICLE%MATERIAL+PARTICLE_HEADER 7   e  H   a   PARTICLE%LAST_MATERIAL+PARTICLE_HEADER 0   ­  H   a   PARTICLE%SQRTKT+PARTICLE_HEADER 5   õ  H   a   PARTICLE%LAST_SQRTKT+PARTICLE_HEADER 5   =  H   a   PARTICLE%N_COLLISION+PARTICLE_HEADER 5     ¤   a   PARTICLE%WRITE_TRACK+PARTICLE_HEADER 1   )  ¤   a   PARTICLE%YES_SAB+PARTICLE_HEADER /   Í  S   a   PARTICLE%CLEAR+PARTICLE_HEADER &      R      CLEAR+PARTICLE_HEADER +   r  V   a   CLEAR%THIS+PARTICLE_HEADER 4   È  X   a   PARTICLE%INITIALIZE+PARTICLE_HEADER +      R      INITIALIZE+PARTICLE_HEADER 0   r  V   a   INITIALIZE%THIS+PARTICLE_HEADER -   È  Z   a   PARTICLE%SET+PARTICLE_HEADER -   "  ^      SET_PARTICLE+PARTICLE_HEADER 2     V   a   SET_PARTICLE%THIS+PARTICLE_HEADER 4   Ö  R   a   SET_PARTICLE%SOURCE+PARTICLE_HEADER 3   (         FOR_ARRAY_DESCRIPTOR+ISO_C_BINDING 8   È  H   a   FOR_ARRAY_DESCRIPTOR%BASE+ISO_C_BINDING 7     H   a   FOR_ARRAY_DESCRIPTOR%LEN+ISO_C_BINDING :   X  H   a   FOR_ARRAY_DESCRIPTOR%OFFSET+ISO_C_BINDING 9      H   a   FOR_ARRAY_DESCRIPTOR%FLAGS+ISO_C_BINDING 8   è  H   a   FOR_ARRAY_DESCRIPTOR%RANK+ISO_C_BINDING =   0  H   a   FOR_ARRAY_DESCRIPTOR%RESERVED1+ISO_C_BINDING ;   x  ²   a   FOR_ARRAY_DESCRIPTOR%DIMINFO+ISO_C_BINDING /   *  v      FOR_DESC_TRIPLET+ISO_C_BINDING 6      H   a   FOR_DESC_TRIPLET%EXTENT+ISO_C_BINDING 4   è  H   a   FOR_DESC_TRIPLET%MULT+ISO_C_BINDING :   0   H   a   FOR_DESC_TRIPLET%LOWERBOUND+ISO_C_BINDING '   x   P       #UNLPOLY+ISO_C_BINDING !   È   y      BANK+BANK_HEADER %   A!  H   a   BANK%WGT+BANK_HEADER %   !     a   BANK%XYZ+BANK_HEADER %   %"     a   BANK%UVW+BANK_HEADER #   Á"  H   a   BANK%E+BANK_HEADER #   	#  H   a   BANK%G+BANK_HEADER     Q#  ®       MG_XS+XS_HEADER '   ÿ#  P   a   MG_XS%MAT_ID+XS_HEADER '   O$     a   MG_XS%SIG_TR+XS_HEADER (   ã$     a   MG_XS%SIG_ABS+XS_HEADER (   w%     a   MG_XS%SIG_CAP+XS_HEADER (   &     a   MG_XS%SIG_FIS+XS_HEADER #   &     a   MG_XS%NU+XS_HEADER $   3'     a   MG_XS%CHI+XS_HEADER )   Ç'  ¬   a   MG_XS%SIG_SCAT+XS_HEADER    s(  @       KEFF+VARIABLES     ³(  @       K_COL+VARIABLES    ó(  @       K_TL+VARIABLES )   3)         FISSION_BANK+BANK_HEADER (   É)         THREAD_BANK+BANK_HEADER %   g*  @       BANK_IDX+BANK_HEADER     §*         XS_MG+XS_HEADER    >+         RANG+RANDOMS !   ×+  Û       RAND_VEC+RANDOMS &   ²,  =      RAND_VEC%SQRT+RANDOMS %   ï,  <      RAND_VEC%COS+RANDOMS %   +-  <      RAND_VEC%SIN+RANDOMS "   g-  t       WGT_MIN+CONSTANTS    Û-  q       BASE+CONSTANTS    L.         COLLISION_MG "   Þ.  =      COLLISION_MG%SIZE !   /  <      COLLISION_MG%INT !   W/  <      COLLISION_MG%SUM    /  V   a   COLLISION_MG%P     é/  ¬       COLLISION_MG_DT %   0  =      COLLISION_MG_DT%SIZE $   Ò0  <      COLLISION_MG_DT%INT $   1  <      COLLISION_MG_DT%SUM "   J1  V   a   COLLISION_MG_DT%P ,    1  @   a   COLLISION_MG_DT%MACRO_MAJOR 