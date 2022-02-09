PRO Transparent_Polygon

   WINDOWSIZE=400

   ; Create some data.
   signal = LoadData(1) - 15
   time = Findgen(N_Elements(signal)) * 6.0 / N_Elements(signal)
   
   ; Create some windows.
   Window, Title='Data Window', XSIZE=WINDOWSIZE, YSIZE=WINDOWSIZE, $
        /FREE, XPOS=0, YPOS=0
   dataWin = !D.Window
   Window, Title='Polygon Window', XSIZE=WINDOWSIZE, YSIZE=WINDOWSIZE, $
        /FREE, XPOS=WINDOWSIZE+10, YPOS=0
   polyWin = !D.Window
   
   ; Draw plot in data window.
   WSet, dataWin
   Plot, time, signal, BACKGROUND=FSC_Color('ivory'), $
      COLOR=FSC_Color('navy'), $
      /NODATA, XTitle='Time', YTitle='Signal Strength'
   OPLOT, time, signal, THICK=2, COLOR=FSC_Color('cornflower blue')
   OPLOT, time, signal, PSYM=2, COLOR=FSC_Color('olive')
   
   ; Take a snapshot.
   background = TVREAD(TRUE=3)
   
   ; Copy data window and draw a polygon in the polygon window.
   WSet, polyWin
   DEVICE, COPY=[0, 0, WINDOWSIZE, WINDOWSIZE, 0, 0, dataWin]
   POLYFILL, [ 0,  6, 6, 0,  0], $
             [-5, -5, 5, 5, -5], /DATA, $
             COLOR=FSC_COLOR('deep pink')
             
   ; Take a snapshot of this window, then delete it.
   foreground = TVREAD(TRUE=3)   
   
   ; Use a 25% transparent alpha.
   alpha = 0.25
   Window, Title='Transparent Polygon Window', XSIZE=WINDOWSIZE, YSIZE=WINDOWSIZE, $
        XPOS=2*WINDOWSIZE+20, YPOS=0, /FREE
   TV, (foreground * alpha) + (1 - alpha) * background, TRUE=3
   
END
