	  �4  z   k820309    h          14.0        ���]                                                                                                           
       tally.f90 TALLY                                                     
                        �                                 
                            @                              
       PARTICLE                   @                               '�                   #N_COORD    #CELL_INSTANCE    #COORD    #LAST_N_COORD    #LAST_CELL    #E    #LAST_E    #G    #LAST_G    #WGT    #MU    #ALIVE    #LAST_UVW    #LAST_WGT    #MATERIAL     #LAST_MATERIAL !   #SQRTKT "   #LAST_SQRTKT #   #N_COLLISION $   #N_CROSS %   #YES_SAB &   #VRC_TRACED '   #CLEAR (   #INITIALIZE +   #SET .                � $                                                              � $                                                             � $                                   
              P             #LOCALCOORD    p          p 
           p 
                                        @  @                              'P              
      #CELL 	   #UNIVERSE 
   #LATTICE    #LATTICE_X    #LATTICE_Y    #LATTICE_Z    #XYZ    #UVW    #DIST    #RESET                �                               	                                                                                                 0                �                               
                                                                                                0                �                                                                                                                               0                �                                                                                                                               0                �                                                                                                                               0                �                                                                                                                               0                 �                                                           
  p          p            p                                       �                                          0                 
  p          p            p                                       �                                   H       	   
   1         �   �                       �                   
     #RESET_COORD    #         H     @                                               #THIS              
                                     P               #LOCALCOORD                 � $                                   (                         � $                                   
       ,                  p          p 
           p 
                                      � $                                  X         
                � $                                  `         
                � $                                   h                         � $                                   l      	                   � $                                  p      
   
                � $                                  x         
                � $                                   �                         � $                                         �                
  p          p            p                                       � $                                  �         
                � $                                    �                         � $                              !     �                         � $                             "     �         
                � $                             #     �         
                � $                              $     �                         � $                              %     �                        � $                              &     �                                                                                                           � $                              '     �                                                                                               1         �   � $                      �      (                  #CLEAR )   #         @     @                            )                    #THIS *                                             *     �              #PARTICLE    1         �   � $                      �      +                  #INITIALIZE ,   #         @     @                            ,                    #THIS -                                             -     �              #PARTICLE    1         �   � $                      �      .                  #SET_PARTICLE /   #         @     @                            /                    #THIS 0   #SOURCE 1                                             0     �              #PARTICLE                                               1     H               #BANK 2                  @  @                          3     '�                    #BASE 4   #LEN 5   #OFFSET 6   #FLAGS 7   #RANK 8   #RESERVED1 9   #DIMINFO :                �                              4                                �                              5                               �                              6                               �                              7                               �                              8                                �                              9     (                          �                               :            0                    #FOR_DESC_TRIPLET ;   p          p            p                                         @  @                         ;     '                    #EXTENT <   #MULT =   #LOWERBOUND >                �                              <                                �                              =                               �                              >                                    @                          ?     '                                       @  @                          2     'H                    #WGT @   #XYZ A   #UVW B   #E C   #G D                �                              @                
                �                              A                             
  p          p            p                                       �                              B                              
  p          p            p                                       �                              C     8          
                �                               D     @                                                          E                                                       0                                             F                                       
               10%         @                                 G                           #THIS H                                             H     X              #CELL I   #         @                                   J                    #THIS K                                             K     �               #UNIVERSE L                     @                         M     '                    #CELL N   #UNIVERSE O   #LATTICE P   #LATTICE_X Q   #LATTICE_Y R   #LATTICE_Z S   #RESET T               �                               N                                                                                                 0                �                               O                                                                                                0                �                               P                                                                                                0                �                               Q                                                                                                0                �                               R                                                                                                0                �                               S                                                                                                0    1         �   �                       �      T                  #RESET_COORD U   #         H                                  U                    #THIS V             
D                                V                    #LOCALCOORD M                     @                          W     '                   #N_COORD X   #COORD Y   #VOL Z   #FLAG [                �                               X                                �                               Y     
                           #LOCALCOORD M   p          p 
           p 
                                      �                              Z     �          
                �                               [                            `                               \                                   &                                           #COORDSTRUCT W                                             ]                   
                &                                                                                     ^                   
                &                                                                                     _                   
                &                                                                                     `                   
                &                                           (         `                                 a                                      #FINDTALLYBIN%SUM b   #FINDTALLYBIN%ABS c   #FINDTALLYBIN%SIZE d   #P e   p          p            p                                                                     b     SUM                                            c     ABS                                            d     SIZE           
                                  e     �             #PARTICLE                      @               @           L     '�                    #UNIV_TYPE f   #UNIV_ID g   #XYZ h   #NCELL i   #R j   #CELL k   #INITIALIZE l                �                               f                                �                               g                               �                              h                             
  p          p            p                                       �                               i                              �                              j            (                 
            &                                                      �                               k            p                             &                                           1         �   �                       �      l                  #INITIALIZE J                     @               @           I     'X                   #CELL_ID m   #IDX n   #UNIV_ID o   #MAT_IDX p   #FILL q   #NSURF r   #FILLTYPE s   #LIST_OF_SURFACE_IDS t   #NEG_SURF_IDX u   #POS_SURF_IDX v   #TRANSLATION w   #OPERAND_FLAG x   #FILL_TYPE y                �                              m                                        �                               n                               �                               o                               �                               p                               �                               q                                �                               r     $                          �                               s     (                        �                               t            0                             &                                                     �                               u            x              	               &                                                     �                               v            �              
               &                                                      �                              w                            
            &                                                        �                               x     P            1         �   �                      �      y                  #CELL_FILL G      �         fn#fn     �   @   J   GEOMETRY_HEADER    �   @   J   CONSTANTS     8  I   J  PARTICLE_HEADER )   �  �      PARTICLE+PARTICLE_HEADER 1     H   a   PARTICLE%N_COORD+PARTICLE_HEADER 7   `  H   a   PARTICLE%CELL_INSTANCE+PARTICLE_HEADER /   �  �   a   PARTICLE%COORD+PARTICLE_HEADER +   T  �      LOCALCOORD+PARTICLE_HEADER 0     �   a   LOCALCOORD%CELL+PARTICLE_HEADER 4   �  �   a   LOCALCOORD%UNIVERSE+PARTICLE_HEADER 3   g  �   a   LOCALCOORD%LATTICE+PARTICLE_HEADER 5     �   a   LOCALCOORD%LATTICE_X+PARTICLE_HEADER 5   �  �   a   LOCALCOORD%LATTICE_Y+PARTICLE_HEADER 5   V  �   a   LOCALCOORD%LATTICE_Z+PARTICLE_HEADER /   �  �   a   LOCALCOORD%XYZ+PARTICLE_HEADER /   �	  �   a   LOCALCOORD%UVW+PARTICLE_HEADER 0   3
  H   a   LOCALCOORD%DIST+PARTICLE_HEADER 1   {
  Y   a   LOCALCOORD%RESET+PARTICLE_HEADER ,   �
  R      RESET_COORD+PARTICLE_HEADER 1   &  X   a   RESET_COORD%THIS+PARTICLE_HEADER 6   ~  H   a   PARTICLE%LAST_N_COORD+PARTICLE_HEADER 3   �  �   a   PARTICLE%LAST_CELL+PARTICLE_HEADER +   b  H   a   PARTICLE%E+PARTICLE_HEADER 0   �  H   a   PARTICLE%LAST_E+PARTICLE_HEADER +   �  H   a   PARTICLE%G+PARTICLE_HEADER 0   :  H   a   PARTICLE%LAST_G+PARTICLE_HEADER -   �  H   a   PARTICLE%WGT+PARTICLE_HEADER ,   �  H   a   PARTICLE%MU+PARTICLE_HEADER /     H   a   PARTICLE%ALIVE+PARTICLE_HEADER 2   Z  �   a   PARTICLE%LAST_UVW+PARTICLE_HEADER 2   �  H   a   PARTICLE%LAST_WGT+PARTICLE_HEADER 2   >  H   a   PARTICLE%MATERIAL+PARTICLE_HEADER 7   �  H   a   PARTICLE%LAST_MATERIAL+PARTICLE_HEADER 0   �  H   a   PARTICLE%SQRTKT+PARTICLE_HEADER 5     H   a   PARTICLE%LAST_SQRTKT+PARTICLE_HEADER 5   ^  H   a   PARTICLE%N_COLLISION+PARTICLE_HEADER 1   �  H   a   PARTICLE%N_CROSS+PARTICLE_HEADER 1   �  �   a   PARTICLE%YES_SAB+PARTICLE_HEADER 4   �  �   a   PARTICLE%VRC_TRACED+PARTICLE_HEADER /   6  S   a   PARTICLE%CLEAR+PARTICLE_HEADER &   �  R      CLEAR+PARTICLE_HEADER +   �  V   a   CLEAR%THIS+PARTICLE_HEADER 4   1  X   a   PARTICLE%INITIALIZE+PARTICLE_HEADER +   �  R      INITIALIZE+PARTICLE_HEADER 0   �  V   a   INITIALIZE%THIS+PARTICLE_HEADER -   1  Z   a   PARTICLE%SET+PARTICLE_HEADER -   �  ^      SET_PARTICLE+PARTICLE_HEADER 2   �  V   a   SET_PARTICLE%THIS+PARTICLE_HEADER 4   ?  R   a   SET_PARTICLE%SOURCE+PARTICLE_HEADER 3   �  �      FOR_ARRAY_DESCRIPTOR+ISO_C_BINDING 8   1  H   a   FOR_ARRAY_DESCRIPTOR%BASE+ISO_C_BINDING 7   y  H   a   FOR_ARRAY_DESCRIPTOR%LEN+ISO_C_BINDING :   �  H   a   FOR_ARRAY_DESCRIPTOR%OFFSET+ISO_C_BINDING 9   	  H   a   FOR_ARRAY_DESCRIPTOR%FLAGS+ISO_C_BINDING 8   Q  H   a   FOR_ARRAY_DESCRIPTOR%RANK+ISO_C_BINDING =   �  H   a   FOR_ARRAY_DESCRIPTOR%RESERVED1+ISO_C_BINDING ;   �  �   a   FOR_ARRAY_DESCRIPTOR%DIMINFO+ISO_C_BINDING /   �  v      FOR_DESC_TRIPLET+ISO_C_BINDING 6   	  H   a   FOR_DESC_TRIPLET%EXTENT+ISO_C_BINDING 4   Q  H   a   FOR_DESC_TRIPLET%MULT+ISO_C_BINDING :   �  H   a   FOR_DESC_TRIPLET%LOWERBOUND+ISO_C_BINDING '   �  P       #UNLPOLY+ISO_C_BINDING !   1  y      BANK+BANK_HEADER %   �  H   a   BANK%WGT+BANK_HEADER %   �  �   a   BANK%XYZ+BANK_HEADER %   �  �   a   BANK%UVW+BANK_HEADER #   *  H   a   BANK%E+BANK_HEADER #   r  H   a   BANK%G+BANK_HEADER    �  q       NONE+CONSTANTS $   +  r       MAX_COORD+CONSTANTS *   �  Z       CELL_FILL+GEOMETRY_HEADER /   �  R   a   CELL_FILL%THIS+GEOMETRY_HEADER +   I  R       INITIALIZE+GEOMETRY_HEADER 0   �  V   a   INITIALIZE%THIS+GEOMETRY_HEADER    �  �       LOCALCOORD     �  �   a   LOCALCOORD%CELL $   C   �   a   LOCALCOORD%UNIVERSE #   �   �   a   LOCALCOORD%LATTICE %   �!  �   a   LOCALCOORD%LATTICE_X %   2"  �   a   LOCALCOORD%LATTICE_Y %   �"  �   a   LOCALCOORD%LATTICE_Z !   |#  Y   a   LOCALCOORD%RESET    �#  R       RESET_COORD !   '$  X   a   RESET_COORD%THIS    $  {       COORDSTRUCT $   �$  H   a   COORDSTRUCT%N_COORD "   B%  �   a   COORDSTRUCT%COORD     �%  H   a   COORDSTRUCT%VOL !   6&  H   a   COORDSTRUCT%FLAG    ~&  �       TALLYCOORD    '  �       TALLYFLUX    �'  �       TALLYPOWER    3(  �       TALLY1    �(  �       TALLY2    K)  �       FINDTALLYBIN !   9*  <      FINDTALLYBIN%SUM !   u*  <      FINDTALLYBIN%ABS "   �*  =      FINDTALLYBIN%SIZE    �*  V   a   FINDTALLYBIN%P )   D+  �       UNIVERSE+GEOMETRY_HEADER 3   �+  H   a   UNIVERSE%UNIV_TYPE+GEOMETRY_HEADER 1   -,  H   a   UNIVERSE%UNIV_ID+GEOMETRY_HEADER -   u,  �   a   UNIVERSE%XYZ+GEOMETRY_HEADER /   -  H   a   UNIVERSE%NCELL+GEOMETRY_HEADER +   Y-  �   a   UNIVERSE%R+GEOMETRY_HEADER .   �-  �   a   UNIVERSE%CELL+GEOMETRY_HEADER 4   �.  X   a   UNIVERSE%INITIALIZE+GEOMETRY_HEADER %   �.        CELL+GEOMETRY_HEADER -   �/  P   a   CELL%CELL_ID+GEOMETRY_HEADER )   ;0  H   a   CELL%IDX+GEOMETRY_HEADER -   �0  H   a   CELL%UNIV_ID+GEOMETRY_HEADER -   �0  H   a   CELL%MAT_IDX+GEOMETRY_HEADER *   1  H   a   CELL%FILL+GEOMETRY_HEADER +   [1  H   a   CELL%NSURF+GEOMETRY_HEADER .   �1  H   a   CELL%FILLTYPE+GEOMETRY_HEADER 9   �1  �   a   CELL%LIST_OF_SURFACE_IDS+GEOMETRY_HEADER 2   2  �   a   CELL%NEG_SURF_IDX+GEOMETRY_HEADER 2   3  �   a   CELL%POS_SURF_IDX+GEOMETRY_HEADER 1   �3  �   a   CELL%TRANSLATION+GEOMETRY_HEADER 2   ;4  H   a   CELL%OPERAND_FLAG+GEOMETRY_HEADER /   �4  W   a   CELL%FILL_TYPE+GEOMETRY_HEADER 