	  %  [   k820309    h          14.0        �R5]                                                                                                           
       surface_header.f90 SURFACE_HEADER                                                    
                                                           
                      @  @                               '�                    #BASE    #LEN    #OFFSET    #FLAGS    #RANK    #RESERVED1 	   #DIMINFO 
                �                                                              �                                                             �                                                             �                                                             �                                                              �                              	     (                          �                               
            0                    #FOR_DESC_TRIPLET    p          p            p                                         @  @                              '                    #EXTENT    #MULT    #LOWERBOUND                 �                                                              �                                                             �                                                                  @                               '                                                                                                                                 
                   
                  �������                                                                                                                                                                                                                                                                                                                  @               @                '�                    #SURF_ID    #SURF_TYPE    #BC    #PARMTRS    #NEIGHBOR_POS    #NEIGHBOR_NEG                 �                                                                      �                                                              �                                                              �                                                            
  p          p            p                                     �                                           H                             &                                                      �                                           �                             &                                                                                                  �                        &                                           #SURFACE                                                           �                        &                                           #SURFACE                                                      �       #SURFACE    #         @                                   !                   #READ_SURF%LEN "   #SURFOBJ #   #LINE $                                              "     LEN           D                                 #     �               #SURFACE              D @                              $                     1 %         @                                %                           #SURF_TYPE &             
                                 &                    1 #         @                                   '                   #READ_BC%LEN (   #SURFLIST )   #LINE *                                              (     LEN           D @                               )            �                       &                                           #SURFACE              D @                              *                     1 %         @                                +                          #FIND_SURF_IDX%SIZE ,   #THIS -   #SURF_ID .                                              ,     SIZE           @                               -            �                       &                                           #SURFACE                                              .                     1 %         @                                 /                          #SURF_NEG_OR_POS%SQRT 0   #SURF_NEG_OR_POS%ABS 1   #THIS 2   #XYZ 3                                              0     SQRT                                            1     ABS                                            2     �               #SURFACE                                              3                   
     p          p            p                          %         @                                 4                    
       #A 5   #B 6   #C 7   #D 8   #XYZ 9   #UVW :             
                                 5     
                
                                 6     
                
                                 7     
                
                                 8     
                
                                 9                   
 
   p          p            p                                    
                                 :                   
    p          p            p                          %         @                                ;                    
       #SURF <   #XYZ =   #UVW >                                              <     �               #SURFACE              
                                 =                   
    p          p            p                                    
                                 >                   
    p          p            p                          %         @                                ?                    
       #SURF @   #XYZ A   #UVW B                                              @     �               #SURFACE              
                                 A                   
    p          p            p                                    
                                 B                   
    p          p            p                          %         @                                C                    
       #SURF D   #XYZ E   #UVW F                                              D     �               #SURFACE              
                                 E                   
    p          p            p                                    
                                 F                   
    p          p            p                          %         @                                G                   
       #SURF_CYLZ%SQRT H   #SURF I   #XYZ J   #UVW K                                              H     SQRT                                            I     �               #SURFACE              
                                 J                   
    p          p            p                                    
                                 K                   
    p          p            p                          %         @                                L                   
       #SURF_SQCZ%MINVAL M   #SURF N   #XYZ O   #UVW P                                              M     MINVAL                                            N     �               #SURFACE              
                                 O                   
    p          p            p                                    
                                 P                   
    p          p            p                          %         @                                Q                   
       #SURF_SPH%SQRT R   #SURF S   #XYZ T   #UVW U                                              R     SQRT                                            S     �               #SURFACE              
                                 T                   
    p          p            p                                    
                                 U                   
    p          p            p                          #         @                                   V                    #SURF W   #XYZ X   #UVW Y   #DIST Z             D @                               W     �               #SURFACE              
  @                              X                   
    p          p            p                                    
  @                              Y                   
    p          p            p                                    D                                Z     
          �   *      fn#fn    �   @   J   CONSTANTS    
  @   J   OMP_LIB 3   J  �      FOR_ARRAY_DESCRIPTOR+ISO_C_BINDING 8   �  H   a   FOR_ARRAY_DESCRIPTOR%BASE+ISO_C_BINDING 7   2  H   a   FOR_ARRAY_DESCRIPTOR%LEN+ISO_C_BINDING :   z  H   a   FOR_ARRAY_DESCRIPTOR%OFFSET+ISO_C_BINDING 9   �  H   a   FOR_ARRAY_DESCRIPTOR%FLAGS+ISO_C_BINDING 8   
  H   a   FOR_ARRAY_DESCRIPTOR%RANK+ISO_C_BINDING =   R  H   a   FOR_ARRAY_DESCRIPTOR%RESERVED1+ISO_C_BINDING ;   �  �   a   FOR_ARRAY_DESCRIPTOR%DIMINFO+ISO_C_BINDING /   L  v      FOR_DESC_TRIPLET+ISO_C_BINDING 6   �  H   a   FOR_DESC_TRIPLET%EXTENT+ISO_C_BINDING 4   
  H   a   FOR_DESC_TRIPLET%MULT+ISO_C_BINDING :   R  H   a   FOR_DESC_TRIPLET%LOWERBOUND+ISO_C_BINDING '   �  P       #UNLPOLY+ISO_C_BINDING    �  @       PX+CONSTANTS #   *  p       INFINITY+CONSTANTS    �  @       PY+CONSTANTS    �  @       PZ+CONSTANTS      @       CYLZ+CONSTANTS    Z  @       SQCZ+CONSTANTS    �  @       SPH+CONSTANTS    �  �       SURFACE       P   a   SURFACE%SURF_ID "   �  H   a   SURFACE%SURF_TYPE    	  H   a   SURFACE%BC     _	  �   a   SURFACE%PARMTRS %   �	  �   a   SURFACE%NEIGHBOR_POS %   �
  �   a   SURFACE%NEIGHBOR_NEG    #  �       SURFACES    �  �       SURFACES_TEMP    U  M       SURF_OBJ    �  r       READ_SURF      <      READ_SURF%LEN "   P  U   a   READ_SURF%SURFOBJ    �  L   a   READ_SURF%LINE $   �  _       SURF_TYPE_CONVERTER .   P  L   a   SURF_TYPE_CONVERTER%SURF_TYPE    �  q       READ_BC      <      READ_BC%LEN !   I  �   a   READ_BC%SURFLIST    �  L   a   READ_BC%LINE    .         FIND_SURF_IDX #   �  =      FIND_SURF_IDX%SIZE #   �  �   a   FIND_SURF_IDX%THIS &   �  L   a   FIND_SURF_IDX%SURF_ID     �  �       SURF_NEG_OR_POS %   e  =      SURF_NEG_OR_POS%SQRT $   �  <      SURF_NEG_OR_POS%ABS %   �  U   a   SURF_NEG_OR_POS%THIS $   3  �   a   SURF_NEG_OR_POS%XYZ    �  ~       SURF_GP    E  @   a   SURF_GP%A    �  @   a   SURF_GP%B    �  @   a   SURF_GP%C      @   a   SURF_GP%D    E  �   a   SURF_GP%XYZ    �  �   a   SURF_GP%UVW    m  l       SURF_PX    �  U   a   SURF_PX%SURF    .  �   a   SURF_PX%XYZ    �  �   a   SURF_PX%UVW    V  l       SURF_PY    �  U   a   SURF_PY%SURF      �   a   SURF_PY%XYZ    �  �   a   SURF_PY%UVW    ?  l       SURF_PZ    �  U   a   SURF_PZ%SURF       �   a   SURF_PZ%XYZ    �  �   a   SURF_PZ%UVW    (  �       SURF_CYLZ    �  =      SURF_CYLZ%SQRT    �  U   a   SURF_CYLZ%SURF    :  �   a   SURF_CYLZ%XYZ    �  �   a   SURF_CYLZ%UVW    b  �       SURF_SQCZ !   �  ?      SURF_SQCZ%MINVAL    #  U   a   SURF_SQCZ%SURF    x  �   a   SURF_SQCZ%XYZ       �   a   SURF_SQCZ%UVW    �          SURF_SPH    !  =      SURF_SPH%SQRT    \!  U   a   SURF_SPH%SURF    �!  �   a   SURF_SPH%XYZ    E"  �   a   SURF_SPH%UVW    �"  n       SURF_SELECT !   G#  U   a   SURF_SELECT%SURF     �#  �   a   SURF_SELECT%XYZ     0$  �   a   SURF_SELECT%UVW !   �$  @   a   SURF_SELECT%DIST 