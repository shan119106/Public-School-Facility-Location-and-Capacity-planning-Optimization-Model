\ Model Lpproblem
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  5996.25 Y11 + 6615 Y12 + 13053.75 Y13 + 9585 Y21 + 10320 Y22 + 12857.5 Y23
   + 13160 Y31 + 9735 Y32 + 15675 Y33 + 7222 Y41 + 7062.5 Y42 + 11400 Y43
   + 6868.75 Y51 + 12250 Y52 + 13725 Y53 + 0.25 Yc1 + 0.25 Yc2 + 0.25 Yc3
   + 0.25 Yc4 + 0.25 Yc5 + 0.075 X1 + 0.05 X1 + 0.075 X1 + 0.15 X1
   + 0.175 X1 + 0.175 X1 + 0.15 X1 + 0.125 X1 + 0.2 X2 + 0.15 X2 + 0.05 X2
   + 0.25 X2 + 0.175 X2 + 0.1 X2 + 0.15 X2 + 0.15 X2 + 0.2 X3 + 0.15 X3
   + 0.15 X3 + 0.175 X3 + 0.05 X3 + 0.075 X3 + 0.025 X3 + 0.025 X3 + 0.3 X4
   + 0.25 X4 + 0.275 X4 + 0.2 X4 + 0.1 X4 + 0.2 X4 + 0.125 X4 + 0.15 X4
   + 0.05 X5 + 0.075 X5 + 0.175 X5 + 0.125 X5 + 0.225 X5 + 0.25 X5
   + 0.225 X5 + 0.2 X5 + 0.25 rho + 25 MinusS + 25 MinusS + 25 MinusS
   + 25 MinusS + 25 MinusS + 25 MinusS + 25 MinusS + 25 MinusS + 25 MinusS
   + 25 MinusS + 25 MinusS + 25 MinusS + 25 MinusS + 25 MinusS + 25 MinusS
   + 25 MinusS + 25 MinusS + 25 MinusS + 25 MinusS + 25 MinusS + 25 MinusS
   + 25 MinusS + 25 MinusS + 25 MinusS + 25 MinusS + 25 MinusS + 25 MinusS
   + 25 MinusS + 25 MinusS + 25 MinusS + 25 MinusS + 25 MinusS + 25 MinusS
   + 25 MinusS + 25 MinusS + 25 MinusS + 25 MinusS + 25 MinusS + 25 MinusS
   + 25 MinusS
Subject To
 R0: Y11 + Y12 + Y13 + Yc1 + Ync1 = 1
 R1: Y21 + Y22 + Y23 + Yc2 + Ync2 = 1
 R2: Y31 + Y32 + Y33 + Yc3 + Ync3 = 1
 R3: Y41 + Y42 + Y43 + Yc4 + Ync4 = 1
 R4: Y51 + Y52 + Y53 + Yc5 + Ync5 = 1
 R5: Yc1 = 0
 R6: Yc2 = 0
 R7: X1 + X2 + X3 + X4 + X5 + V1 = 139
 R8: X1 + X2 + X3 + X4 + X5 + V2 = 222
 R9: X1 + X2 + X3 + X4 + X5 + V3 = 133
 R10: X1 + X2 + X3 + X4 + X5 + V4 = 121
 R11: X1 + X2 + X3 + X4 + X5 + V5 = 131
 R12: X1 + X2 + X3 + X4 + X5 + V6 = 169
 R13: X1 + X2 + X3 + X4 + X5 + V7 = 162
 R14: X1 + X2 + X3 + X4 + X5 + V8 = 111
 R15: - 250 Y11 - 280 Y12 - 300 Y13 - 200 Ync1 + X1 + X1 + X1 + X1 + X1
   + X1 + X1 + X1 + U1 = 0
 R16: - 300 Y21 - 320 Y22 - 370 Y23 - 250 Ync2 + X2 + X2 + X2 + X2 + X2
   + X2 + X2 + X2 + U2 = 0
 R17: - 300 Y31 - 330 Y32 - 380 Y33 - 250 Ync3 + X3 + X3 + X3 + X3 + X3
   + X3 + X3 + X3 + U3 = 0
 R18: - 300 Y41 - 320 Y42 - 350 Y43 - 220 Ync4 + X4 + X4 + X4 + X4 + X4
   + X4 + X4 + X4 + U4 = 0
 R19: - 250 Y51 - 300 Y52 - 350 Y53 - 226 Ync5 + X5 + X5 + X5 + X5 + X5
   + X5 + X5 + X5 + U5 = 0
 R20: X1 <= 139
 R21: X1 <= 222
 R22: X1 <= 133
 R23: X1 <= 121
 R24: X1 <= 131
 R25: X1 <= 169
 R26: X1 <= 162
 R27: X1 <= 111
 R28: X2 <= 0
 R29: X2 <= 222
 R30: X2 <= 133
 R31: X2 <= 0
 R32: X2 <= 131
 R33: X2 <= 169
 R34: X2 <= 162
 R35: X2 <= 111
 R36: X3 <= 0
 R37: X3 <= 222
 R38: X3 <= 133
 R39: X3 <= 121
 R40: X3 <= 131
 R41: X3 <= 169
 R42: X3 <= 162
 R43: X3 <= 111
 R44: X4 <= 0
 R45: X4 <= 0
 R46: X4 <= 0
 R47: X4 <= 0
 R48: X4 <= 131
 R49: X4 <= 0
 R50: X4 <= 162
 R51: X4 <= 111
 R52: X5 <= 139
 R53: X5 <= 222
 R54: X5 <= 133
 R55: X5 <= 121
 R56: X5 <= 0
 R57: X5 <= 0
 R58: X5 <= 0
 R59: X5 <= 0
 R60: X1 - plusS + MinusS = 1
 R61: X1 - plusS + MinusS = 1
 R62: X1 - plusS + MinusS = 1
 R63: X1 - plusS + MinusS = 1
 R64: X1 - plusS + MinusS = 1
 R65: X1 - plusS + MinusS = 1
 R66: X1 - plusS + MinusS = 1
 R67: X1 - plusS + MinusS = 1
 R68: X2 - plusS + MinusS = 1
 R69: X2 - plusS + MinusS = 1
 R70: X2 - plusS + MinusS = 1
 R71: X2 - plusS + MinusS = 1
 R72: X2 - plusS + MinusS = 1
 R73: X2 - plusS + MinusS = 1
 R74: X2 - plusS + MinusS = 1
 R75: X2 - plusS + MinusS = 1
 R76: X3 - plusS + MinusS = 1
 R77: X3 - plusS + MinusS = 1
 R78: X3 - plusS + MinusS = 1
 R79: X3 - plusS + MinusS = 1
 R80: X3 - plusS + MinusS = 1
 R81: X3 - plusS + MinusS = 1
 R82: X3 - plusS + MinusS = 1
 R83: X3 - plusS + MinusS = 1
 R84: X4 - plusS + MinusS = 1
 R85: X4 - plusS + MinusS = 1
 R86: X4 - plusS + MinusS = 1
 R87: X4 - plusS + MinusS = 1
 R88: X4 - plusS + MinusS = 1
 R89: X4 - plusS + MinusS = 1
 R90: X4 - plusS + MinusS = 1
 R91: X4 - plusS + MinusS = 1
 R92: X5 - plusS + MinusS = 1
 R93: X5 - plusS + MinusS = 1
 R94: X5 - plusS + MinusS = 1
 R95: X5 - plusS + MinusS = 1
 R96: X5 - plusS + MinusS = 1
 R97: X5 - plusS + MinusS = 1
 R98: X5 - plusS + MinusS = 1
 R99: X5 - plusS + MinusS = 1
 R100: - 50 Y11 - 56 Y12 - 60 Y13 - 40 Ync1 + U1 <= 0
 R101: - 60 Y21 - 64 Y22 - 74 Y23 - 50 Ync2 + U2 <= 0
 R102: - 60 Y31 - 66 Y32 - 76 Y33 - 50 Ync3 + U3 <= 0
 R103: - 60 Y41 - 64 Y42 - 70 Y43 - 44 Ync4 + U4 <= 0
 R104: - 50 Y51 - 60 Y52 - 70 Y53 - 45.2 Ync5 + U5 <= 0
 R105: V1 <= 41.7
 R106: V2 <= 66.6
 R107: V3 <= 39.9
 R108: V4 <= 36.3
 R109: V5 <= 39.3
 R110: V6 <= 50.7
 R111: V7 <= 48.6
 R112: V8 <= 33.3
 R113: 0.3 X1 + 0.8 X2 + 0.8 X3 + 1.2 X4 + 0.2 X5 - 139 rho <= 0
 R114: 0.2 X1 + 0.6 X2 + 0.6 X3 + X4 + 0.3 X5 - 222 rho <= 0
 R115: 0.3 X1 + 0.2 X2 + 0.6 X3 + 1.1 X4 + 0.7 X5 - 133 rho <= 0
 R116: 0.6 X1 + X2 + 0.7 X3 + 0.8 X4 + 0.5 X5 - 121 rho <= 0
 R117: 0.7 X1 + 0.7 X2 + 0.2 X3 + 0.4 X4 + 0.9 X5 - 131 rho <= 0
 R118: 0.7 X1 + 0.4 X2 + 0.3 X3 + 0.8 X4 + X5 - 169 rho <= 0
 R119: 0.6 X1 + 0.6 X2 + 0.1 X3 + 0.5 X4 + 0.9 X5 - 162 rho <= 0
 R120: 0.5 X1 + 0.6 X2 + 0.1 X3 + 0.6 X4 + 0.8 X5 - 111 rho <= 0
Bounds
Binaries
 Y11 Y12 Y13 Y21 Y22 Y23 Y31 Y32 Y33 Y41 Y42 Y43 Y51 Y52 Y53 Yc1 Yc2 Yc3
 Yc4 Yc5 Ync1 Ync2 Ync3 Ync4 Ync5
End
