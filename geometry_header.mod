	  ×<     k820309    h          14.0        ©vÊ\                                                                                                           
       geometry_header.f90 GEOMETRY_HEADER                                                    
                                                           
                                                          
                         @               @                'Ø                    #SURF_ID    #SURF_TYPE    #BC    #PARMTRS    #NEIGHBOR_POS 	   #NEIGHBOR_NEG 
                                                                                                                                                                                                                                                                              
  p          p            p                                                                    	            H                             &                                                                                     
                                         &                                                          @  @                               'Ø                    #BASE    #LEN    #OFFSET    #FLAGS    #RANK    #RESERVED1    #DIMINFO                                                                                                                                                                                                                                                                                                                                                                       (                                                                     0                    #FOR_DESC_TRIPLET    p          p            p                                         @  @                              '                    #EXTENT    #MULT    #LOWERBOUND                                                                                                                                                                                                              @                               '                                                                                                                           1                                                                                                   2                                                                                                   3                                                 
                   
                  ÿÿÿÿÿÿï        #         @                                                      #SURF    #XYZ    #UVW    #DIST                                                     Ø               #SURFACE              
                                                    
    p          p            p                                    
                                                    
    p          p            p                                                                          
       %         @                                !                          #SURF_NEG_OR_POS%ABS "   #SURF_NEG_OR_POS%SQRT #   #THIS $   #XYZ %                                               "     ABS                                             #     SQRT                                            $     Ø               #SURFACE                                              %                   
     p          p            p                                   @ `                               &            Ø                        &                                           #SURFACE                                                '                                                        (                                                        )                                                        *                                                        +                              @               @           ,     'X                   #CELL_ID -   #IDX .   #UNIV_ID /   #MAT_IDX 0   #FILL 1   #NSURF 2   #FILLTYPE 3   #LIST_OF_SURFACE_IDS 4   #NEG_SURF_IDX 5   #POS_SURF_IDX 6   #TRANSLATION 7   #OPERAND_FLAG 8   #FILL_TYPE 9                                              -                                                                       .                                                              /                                                              0                                                              1                                                               2     $                                                         3     (                                                       4            0                             &                                                                                    5            x              	               &                                                                                    6            À              
               &                                                                                    7                            
            &                                                                                       8     P            1         À                               9                  #CELL_FILL :   %         @                                 :                           #THIS ;             D @                              ;     X              #CELL ,                                              <            X                       &                                           #CELL ,                                              =            X                       &                                           #CELL ,                                               >     X      #CELL ,                     @               @           ?     '¸                    #UNIV_TYPE @   #UNIV_ID A   #XYZ B   #NCELL C   #R D   #CELL E   #INITIALIZE F                                               @                                                               A                                                             B                             
  p          p            p                                                                      C                                                            D            (                 
            &                                                                                     E            p                             &                                           1         À                                F                  #INITIALIZE G   #         @                                   G                    #THIS H             D                                H     ¸               #UNIVERSE ?            @ `                               I            ¸                        &                                           #UNIVERSE ?                                              J            ¸                        &                                           #UNIVERSE ?                                               K     ¸       #UNIVERSE ?                                              L     ¸       #UNIVERSE ?                     @               @           M     'À                    #LAT_ID N   #LAT_TYPE O   #XYZ P   #N_XYZ Q   #LAT R   #PITCH S                                               N                                                               O                                                             P                             
  p          p            p                                                                      Q                                p          p            p                                                                    R            0                             &                   &                   &                                                                                      S            ¨                 
  p          p            p                                   @ `                               T            À                        &                                           #LATTICE M                                              U            À                        &                                           #LATTICE M                                               V     À       #LATTICE M                                              W     À       #LATTICE M                                             X                   
                &                                           %         @                                 Y                          #FIND_UNIV_IDX%SIZE Z   #THIS [   #UNIV_ID \                                              Z     SIZE           @                               [            ¸                       & p                                           #UNIVERSE ?                                              \            %         @                                 ]                          #FIND_CELL_IDX%SIZE ^   #THIS _   #CELL_ID `                                              ^     SIZE           @                               _            X                      &                                           #CELL ,                                             `                     1 %         @                                 a                          #FIND_LAT_IDX%SIZE b   #THIS c   #LAT_ID d                                              b     SIZE           @                               c            À                       &                                           #LATTICE M                                              d            (         `                                 e                                      #LATTICE_COORD%CEILING f   #LATTICE_COORD%FLOOR g   #LATTICE_COORD%MOD h   #LATTICE_COORD%SIZE i   #LATTICE_COORD%REAL j   #THIS k   #XYZ0 l   p          p            p                                                                     f     CEILING                                            g     FLOOR                                            h     MOD                                            i     SIZE                             @              j     REAL                                            k     À               #LATTICE M                                             l                   
     p          p            p                          %         @                                m                          #IN_THE_LIST_UNIV%SIZE n   #THIS o   #NAME p                                              n     SIZE           @                              o            ¸                       &                                           #UNIVERSE ?                                              p            %         @                                q                          #IN_THE_LIST_LAT%SIZE r   #THIS s   #NAME t                                              r     SIZE           @                              s            À                       &                                           #LATTICE M                                              t            (         `                                 u                                  
    #GET_LOCAL_XYZ%CEILING v   #GET_LOCAL_XYZ%MOD w   #GET_LOCAL_XYZ%REAL x   #THIS y   #UPPER_XYZ z   #I_XYZ {   p          p            p                                                                     v     CEILING                                            w     MOD                             @              x     REAL                                            y     À               #LATTICE M                                             z                   
     p          p            p                                                                     {                        p          p            p                          #         @                                  |                   #CELL_DISTANCE%SIZE }   #THIS ~   #XYZ    #UVW    #SURFLIST    #D_SURF    #IDX_SURF                                               }     SIZE                                            ~     X              #CELL ,              @                                                 
      p          p            p                                     @                                                 
 !    p          p            p                                    D @                                           Ø                       &                                           #SURFACE              D                                     
                 D                                             #         @                                                       #THIS    #SURFLIST    #XYZ    #UVW    #I_XYZ    #D_SURF    #IDX_SURF              
                                       À              #LATTICE M             
@ @                                           Ø       "               &                                           #SURFACE              
                                                    
 $   p          p            p                                    
@ @                                                 
 %   p          p            p                                    
                                                      #   p          p            p                                    D                                     
                 D                                             %         @                                                          #CELL_CONTAINS_XYZ%SIZE    #C    #XYZ                                                    SIZE           
                                       X             #CELL ,             
@ @                                                 
 '   p          p            p                                 ,      fn#fn    Ì   @   J   CONSTANTS      @   J   SURFACE_HEADER    L  @   J   OMP_LIB '     ¥       SURFACE+SURFACE_HEADER /   1  P   a   SURFACE%SURF_ID+SURFACE_HEADER 1     H   a   SURFACE%SURF_TYPE+SURFACE_HEADER *   É  H   a   SURFACE%BC+SURFACE_HEADER /        a   SURFACE%PARMTRS+SURFACE_HEADER 4   ­     a   SURFACE%NEIGHBOR_POS+SURFACE_HEADER 4   A     a   SURFACE%NEIGHBOR_NEG+SURFACE_HEADER 3   Õ         FOR_ARRAY_DESCRIPTOR+ISO_C_BINDING 8   u  H   a   FOR_ARRAY_DESCRIPTOR%BASE+ISO_C_BINDING 7   ½  H   a   FOR_ARRAY_DESCRIPTOR%LEN+ISO_C_BINDING :     H   a   FOR_ARRAY_DESCRIPTOR%OFFSET+ISO_C_BINDING 9   M  H   a   FOR_ARRAY_DESCRIPTOR%FLAGS+ISO_C_BINDING 8     H   a   FOR_ARRAY_DESCRIPTOR%RANK+ISO_C_BINDING =   Ý  H   a   FOR_ARRAY_DESCRIPTOR%RESERVED1+ISO_C_BINDING ;   %  ²   a   FOR_ARRAY_DESCRIPTOR%DIMINFO+ISO_C_BINDING /   ×  v      FOR_DESC_TRIPLET+ISO_C_BINDING 6   M  H   a   FOR_DESC_TRIPLET%EXTENT+ISO_C_BINDING 4     H   a   FOR_DESC_TRIPLET%MULT+ISO_C_BINDING :   Ý  H   a   FOR_DESC_TRIPLET%LOWERBOUND+ISO_C_BINDING '   %	  P       #UNLPOLY+ISO_C_BINDING (   u	  q       FILL_MATERIAL+CONSTANTS (   æ	  q       FILL_UNIVERSE+CONSTANTS '   W
  q       FILL_LATTICE+CONSTANTS #   È
  p       INFINITY+CONSTANTS +   8  n       SURF_SELECT+SURFACE_HEADER 0   ¦  U   a   SURF_SELECT%SURF+SURFACE_HEADER /   û     a   SURF_SELECT%XYZ+SURFACE_HEADER /        a   SURF_SELECT%UVW+SURFACE_HEADER 0   #  @   a   SURF_SELECT%DIST+SURFACE_HEADER /   c         SURF_NEG_OR_POS+SURFACE_HEADER 3   ù  <      SURF_NEG_OR_POS%ABS+SURFACE_HEADER 4   5  =      SURF_NEG_OR_POS%SQRT+SURFACE_HEADER 4   r  U   a   SURF_NEG_OR_POS%THIS+SURFACE_HEADER 3   Ç     a   SURF_NEG_OR_POS%XYZ+SURFACE_HEADER (   [         SURFACES+SURFACE_HEADER    ô  @       N_UNIVERSE    4  @       N_CELL    t  @       N_SURFACE    ´  @       BASE_UNIV    ô  @       ISIZE    4        CELL    F  P   a   CELL%CELL_ID      H   a   CELL%IDX    Þ  H   a   CELL%UNIV_ID    &  H   a   CELL%MAT_IDX    n  H   a   CELL%FILL    ¶  H   a   CELL%NSURF    þ  H   a   CELL%FILLTYPE )   F     a   CELL%LIST_OF_SURFACE_IDS "   Ú     a   CELL%NEG_SURF_IDX "   n     a   CELL%POS_SURF_IDX !        a   CELL%TRANSLATION "     H   a   CELL%OPERAND_FLAG    Þ  W   a   CELL%FILL_TYPE    5  Z       CELL_FILL      R   a   CELL_FILL%THIS    á         CELLS    w         CELLS_TEMP      J       OBJ_CELL    W  ¡       UNIVERSE #   ø  H   a   UNIVERSE%UNIV_TYPE !   @  H   a   UNIVERSE%UNIV_ID         a   UNIVERSE%XYZ    $  H   a   UNIVERSE%NCELL    l     a   UNIVERSE%R          a   UNIVERSE%CELL $     X   a   UNIVERSE%INITIALIZE    ì  R       INITIALIZE     >  V   a   INITIALIZE%THIS             UNIVERSES    .         UNIVERSES_TEMP    È  N       OBJ_UNIV      N       UNIVPTR    d         LATTICE    ö  H   a   LATTICE%LAT_ID !   >   H   a   LATTICE%LAT_TYPE          a   LATTICE%XYZ    "!     a   LATTICE%N_XYZ    ¾!  Ä   a   LATTICE%LAT    "     a   LATTICE%PITCH    #         LATTICES    ·#         LATTICES_TEMP    P$  M       OBJ_LAT    $  M       LAT_PTR    ê$         SGRID    v%         FIND_UNIV_IDX #   õ%  =      FIND_UNIV_IDX%SIZE #   2&     a   FIND_UNIV_IDX%THIS &   Ð&  @   a   FIND_UNIV_IDX%UNIV_ID    '         FIND_CELL_IDX #   '  =      FIND_CELL_IDX%SIZE #   Ì'     a   FIND_CELL_IDX%THIS &   b(  L   a   FIND_CELL_IDX%CELL_ID    ®(  }       FIND_LAT_IDX "   +)  =      FIND_LAT_IDX%SIZE "   h)     a   FIND_LAT_IDX%THIS $   *  @   a   FIND_LAT_IDX%LAT_ID    A*  3      LATTICE_COORD &   t+  @      LATTICE_COORD%CEILING $   ´+  >      LATTICE_COORD%FLOOR "   ò+  <      LATTICE_COORD%MOD #   .,  =      LATTICE_COORD%SIZE #   k,  =      LATTICE_COORD%REAL #   ¨,  U   a   LATTICE_COORD%THIS #   ý,     a   LATTICE_COORD%XYZ0 !   -         IN_THE_LIST_UNIV &   .  =      IN_THE_LIST_UNIV%SIZE &   M.     a   IN_THE_LIST_UNIV%THIS &   ç.  @   a   IN_THE_LIST_UNIV%NAME     '/  ~       IN_THE_LIST_LAT %   ¥/  =      IN_THE_LIST_LAT%SIZE %   â/     a   IN_THE_LIST_LAT%THIS %   {0  @   a   IN_THE_LIST_LAT%NAME    »0        GET_LOCAL_XYZ &   Í1  @      GET_LOCAL_XYZ%CEILING "   2  <      GET_LOCAL_XYZ%MOD #   I2  =      GET_LOCAL_XYZ%REAL #   2  U   a   GET_LOCAL_XYZ%THIS (   Û2     a   GET_LOCAL_XYZ%UPPER_XYZ $   o3     a   GET_LOCAL_XYZ%I_XYZ    4  ¤       CELL_DISTANCE #   §4  =      CELL_DISTANCE%SIZE #   ä4  R   a   CELL_DISTANCE%THIS "   65     a   CELL_DISTANCE%XYZ "   Ê5     a   CELL_DISTANCE%UVW '   ^6     a   CELL_DISTANCE%SURFLIST %   ÷6  @   a   CELL_DISTANCE%D_SURF '   77  @   a   CELL_DISTANCE%IDX_SURF    w7         LAT_DISTANCE "   8  U   a   LAT_DISTANCE%THIS &   c8     a   LAT_DISTANCE%SURFLIST !   ü8     a   LAT_DISTANCE%XYZ !   9     a   LAT_DISTANCE%UVW #   $:     a   LAT_DISTANCE%I_XYZ $   ¸:  @   a   LAT_DISTANCE%D_SURF &   ø:  @   a   LAT_DISTANCE%IDX_SURF "   8;  |       CELL_CONTAINS_XYZ '   ´;  =      CELL_CONTAINS_XYZ%SIZE $   ñ;  R   a   CELL_CONTAINS_XYZ%C &   C<     a   CELL_CONTAINS_XYZ%XYZ 