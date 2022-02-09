;*****************************************************************************************************
;+
; NAME:
;       IMGWIN
;
; PURPOSE:
;
;       This is an image display routine to allow the user to interact with
;       an image. Moving the cursor in the window gives the value and location
;       inside the image. The image is enclosed in a SCALEIMAGE object, so many of 
;       the input keywords are used to set the scaling parameters for that object. 
;       If axes are requested, an IMGAXIS object is added to the SCALEIMAGE object. 
;       Other keywords are used to set up the axis object. The entire window can
;       be saved to file in various formats. The image can be scaled interactively,
;       image colors can be changed, and various image and axes properties can be
;       accessed directly from the File menu of the image window.
;       
; AUTHORS:
;
;        FANNING SOFTWARE CONSULTING   BURRIDGE COMPUTING
;        1645 Sheely Drive             18 The Green South
;        Fort Collins                  Warborough, Oxon
;        CO 80526 USA                  OX10 7DN, ENGLAND
;        Phone: 970-221-0438           Phone: +44 (0)1865 858279
;        E-mail: davidf@dfanning.com   E-mail: davidb@burridgecomputing.co.uk
;
; CATEGORY:
;
;       Graphics display.
;
; SYNTAX:
;
;       ImgWin, image
;       
; ARGUMENTS:
; 
;     image    An 8-bit or 24-bit image array, or the name of an image file that
;              can be opened with READ_IMAGE. 
;       
; INPUT_KEYWORDS:
;
;     AXES:        Set this keyword to draw a set of axes around the image.
;
;     BACKGROUND:  This keyword specifies the name of a background color. By default, 'ivory'.
;
;     BETA:        The beta factor in a Hyperpolic Sine stretch. Default is 3.0.
;
;     BOTTOM:      The lowest value of the scaled image.
;     
;     BREWER:      Set if the color table index number (CT) is the index of a Brewer color table.
;                  To use Brewer color tables, the file fsc_brewer.tbl must be in your IDL path.
;                  
;     COLOR:       Set this keyword to the name of the axis color. The name is passed to FSC_Color 
;                  for processing. Default is "charcoal". Used only if AXES is set.
;
;     CTINDEX:     The index of the color table to use to display the image. Applies only to 
;                  2D image arrays. By default, 0, gray scale. If set to -1, uses current color table.
;
;     EXPONENT:    The logarithm exponent in a logarithmic stretch. Default is 4.0.
;
;     GAMMA:       The gamma factor in a gamma stretch. Default is 1.5.
;
;     KEEP_ASPECT:  Normally, the image will be resized to fit the specified position in the 
;                   window. If you prefer, you can force the image to maintain its aspect ratio 
;                   in the window (although not its natural size) by setting this keyword.
;                   The image width is fitted first. If, after setting the image width, the image 
;                   height is too big for the window, then the image height is fitted into the window. 
;                   The appropriate values of the POSITION keyword are honored during this fitting 
;                   process. Once a fit is made, the POSITION coordiates are re-calculated to center 
;                   the image in the window. You can recover these new position coordinates as the 
;                   output from the POSITION keyword. This keyword is turned ON by default. In other 
;                   words, to allow free positioning, set KEEP_ASPECT=0. Note that if this keyword is
;                   set, and the XSIZE and YSIZE keywords are undefined, that the window will have the
;                   same aspect ratio as the image.
;                  
;     MEAN:         The mean factor in a logarithmic stretch. Default is 0.5.
;
;     MISSING_COLOR: The name of the missing color. The default is "gray".
;
;     MISSING_VALUE: The number that represents missing value in the image.
;
;     NCOLORS:       The number of colors to scale the data into, as in this: (Default: 256)
;
;                       displayImage = BYTSCL(image, MIN=self.sclmin, MAX=self.sclmax, TOP=self.ncolors-1)
;
;     NOINTERP:      Setting this keyword disables the default bilinear
;                    interpolation done to the image when it is resized. Nearest
;                    neighbor interpolation is done instead. This is preferred
;                    when you do not wish to change the pixel values of the image.
;                    This keyword is turned ON by default. In other words, to allow
;                    interpolation, set NOINTERP=0.
;
;     POSITION:      The position of the image in the display window. The position is given
;                    as a four-element array in normalized (0 to 1) coordinates of the form
;                    [x0, y0, x1, y1], where (x0,y0) is the lower-left corner of the image and
;                    (x1,y1) is the upper-right corner of the image. If the KEEP_ASPECT keyword
;                    is set, the image will be located within the specified POSITION in a way
;                    that preserves the aspect ratio of the image. The default is [0.075, 0.075, 0.925, 0.925].
;                    
;                    
;     SCALETYPE:      The type of scaling performed prior to display. Default is 0, linear scaling.
;
;           Number    Type of Stretch
;             0         Linear         scaled = BytScl(image, MIN=minThresh, MAX=maxThresh)
;             1         Gamma          scaled = GmaScl(image, MIN=minThresh, MAX=maxThresh, Gamma=gamma)
;             2         Log            scaled = LogScl(image, MIN=minThresh, MAX=maxThresh, Mean=mean, Exponent=exponent)
;             3         Asinh          scaled = AsinhScl(image, MIN=minThresh, MAX=maxThresh, Beta=beta)
;             4         Linear 2%      A linear stretch, with 2 percent of pixels clipped at both the top and bottom
;             5         Square Root    A linear stretch of the square root histogram of the image values.
;             6         Equalization   A linear stretch of the histogram equalized image histogram.
;             7         Gaussian       A Gaussian normal function is applied to the image histogram.
;
;     SCLMIN:         The image data is scaled between SCLMIN and SCLMAX before display. Default = 0.
;
;     SCLMAX:         The image data is scaled between SCLMIN and SCLMAX before display. Default = 255.
;
;     SIGMA:          The sigma scale factor for Gaussian scaling. Default is 1.0.
;     
;     XRANGE:         If the AXES keyword is set, this keyword is a two-element vector
;                     giving the X axis range. By default, [0, size of image in X].
;                    
;     XSIZE:          The X size of the initial image window. If undefined, appoximately 600 pixels.
;                     (Acutally size determined by the aspect ratio of the image.)
;                    
;     XTICKFORMAT:    The tick formatting for the X axis, if the AXES keyword is set.
;                    
;     XTILE:          The title of the X axis, if the AXES keyword is set.
;
;     YRANGE:         If the AXES keyword is set, this keyword is a two-element vector
;                     giving the Y axis range. By default, [0, size of image in Y].
;                    
;     YSIZE:          The Y size of the initial image window. If undefined, appoximately 600 pixels.
;                     (Acutally size determined by the aspect ratio of the image.)
;
;     YTICKFORMAT:    The tick formatting for the Y axis, if the AXES keyword is set.
;                    
;     YTILE:          The title of the Y axis, if the AXES keyword is set.
;     
; OUTPUT_KEYWORDS:
; 
;     OUTIMAGE:        The image object that is created inside the program.
;
; MODIFICATION_HISTORY:
;
;       Written by: David W. Fanning, 12 October 2008.
;-
;*****************************************************************************************************
PRO ImgWin::CreateStatusBar

; The purpose of this method is to create the statusbar for the program.

   ; Create a statusbar object.
   self._statusbar = OBJ_NEW('STATUSBAR', self, Name='Statusbar', /Align_Left)


END
;*****************************************************************************************************



PRO ImgWin::EventHandler, event

; This is the main event handler for the text program. All widget objects
; causing events have been named so we can branch on EVENT.NAME.

   ; Branch on the object name. Check the GUI method to see which names apply
   ; to which widgets. The event names are listed in alphabetical order in the
   ; following CASE statement.
   CASE StrUpCase(event.name) OF

      ; Allows the user to set the axes properties interactively.
      'AXES_PROPERTIES': BEGIN
            self.theDrawWidget -> SetWindow
            self.theAxes -> ControlPanel
        END

      ; Exits the program and destroys the TLB self object.
      'EXIT' : OBJ_DESTROY, self
      
      ; Allows the user to change image colors.
      'IMAGE_COLORS': BEGIN
      
            ; Make sure this is the current graphics widnow.
            self.theDrawWidget -> SetWindow
            
            ; Load the current image colors so that XCOLORS starts
            ; up with the proper colors.
            self.theImage -> GetProperty, COLOR_OBJECT=colors
            colors -> GetProperty, COLORPALETTE=palette
            TVLCT, palette
            
            ; Put the XCOLORS palette next to this window.
            self -> GetProperty, XOFFSET=xoffset, YOFFSET=yoffset, XSIZE=xsize, YSIZE=ysize
            colors -> XColors, XOFFSET=xoffset + xsize + 15, YOFFSET=yoffset + 100, GROUP_LEADER=self
                
        END

      ; Allows the user to set the image properties interactively.
      'IMAGE_PROPERTIES': BEGIN
            self.theDrawWidget -> SetWindow
            self.theImage -> ControlPanel
        END

      ; This is where draw widget events are handled.
      'IMGWIN_DRAWWIDGET': BEGIN
 
         ; Make this draw widget the current window.
         event.id -> SetWindow

         ; Get the image associated with this draw widget.
         imageObject = event.ID -> Get(Position=0)

         ; Find the RGB value of the image at this location and display in statusbar.
         value = imageObject -> Pixel_To_Value(event.x, event.y, Inside=inside, XPixel=xpix, YPixel=ypix, $
            XData=xdata, YData=ydata)
         IF inside EQ 0 THEN self._statusbar -> SetProperty, Text='Outside Image' ELSE $
         BEGIN
            imageObject -> GetProperty, N_DIMENSIONS=ndims
            IF ndims EQ 2 THEN BEGIN
                s = 'Value: ' + StrTrim(value,2) 
                s = s + '    Image Coordinate: (' + StrTrim(xpix,2) + ', ' + StrTrim(ypix,2) + ')' + $
                '    Image Location: (' + String(xdata,Format='(F0.2)') + ', ' + String(ydata,Format='(F0.2)') + ')'                
            ENDIF ELSE BEGIN
                s = 'R: ' + StrTrim(value[0],2) + '  G: ' + StrTrim(value[1],2) + '  B: ' + StrTrim(value[2],2)
                s = s + '    Image Coordinate: (' + StrTrim(xpix,2) + ', ' + StrTrim(ypix,2) + ')' + $
                '    Image Location: (' + String(xdata,Format='(F0.2)') + ', ' + String(ydata,Format='(F0.2)') + ')'
            ENDELSE
            self._statusbar -> SetProperty, Text=s
         ENDELSE
         END
                         
       ; This is a resize event from the TLB.
       'IMGWIN_TLB': BEGIN
            
            ; If KEEP_ASPECT is set for the image, then we will constain the draw widget
            ; to have the same aspect ratio.
            self.theImage -> GetProperty, ASPECT_RATIO=aspect, KEEP_ASPECT=keep_aspect
            sizes = Get_Screen_Size()
            IF sizes[0] GT 2000 THEN sizes[0] = sizes[0] / 2 ; Dual monitor problem
            IF keep_aspect THEN BEGIN
                IF aspect GE 1 THEN BEGIN
                    ysize = event.y < (sizes[1] - 95)
                    xsize = ysize / aspect                  
                ENDIF ELSE BEGIN
                    xsize = event.x < (sizes[0] - 8)
                    ysize = xsize * aspect
                ENDELSE
                IF self.full_resolution EQ 0 THEN $
                    self.theDrawWidget -> SetProperty, XSIZE=xsize, YSIZE=ysize
            ENDIF ELSE BEGIN
                IF self.full_resolution EQ 0 THEN $
                    self.theDrawWidget -> SetProperty, XSIZE=event.x < (sizes[0] - 8), YSIZE=event.y < (sizes[1] - 95)
            ENDELSE
            
            ; Redraw the image in the newly-resized draw widget.
            self.theDrawWidget -> Draw
            
            ; Update the size of the statusbar to reflect the size of the draw widget.
            self._statusbar -> Resize, self.theDrawWidget
        END
        
      ; The SAVE AS... buttons come here. Just get the type of file out of the UVALUE of the
      ; button object and tell the draw widget to create a file of this type.
      'SAVE_WINDOW': BEGIN
     
            ; Call the OUTPUT method on the draw widget.
            event.ID -> GetProperty, UVALUE=fileType
            self.theDrawWidget -> Output, TYPE=fileType, FILENAME='imgwin'
            END
            
      ; This allows the user to scale the image interactively. XSTRETCH is called, and stretch
      ; events are sent to the XSTRETCH_NOTIFICATION method.
      'SCALE_IMAGE': BEGIN

            ; Get the current stretch parameters from the image so we can configure
            ; XSTRETCH properly.
            self.theImage -> GetProperty, $
               BETA=beta, $
               BOTTOM=bottom, $
               EXPONENT=exponent, $
               GAMMA=gamma, $
               IMAGE=image, $
               MEAN=mean, $
               NCOLORS=ncolors, $
               NEGATIVE=negative, $
               SCALETYPE=scaletype, $
               SCLMIN=sclmin, $
               SCLMAX=sclmax, $
               SIGMA=sigma

            ; Make sure you are drawing in the right window. Start XSTRETCH.
            self.theDrawWidget -> SetWindow
            XStretch, image, GROUP_LEADER=self->GetID(), /NO_WINDOW, $
               BETA=beta, $
               BOTTOM=bottom, $
               EXPONENT=exponent, $
               GAMMA=gamma, $
               MEAN=mean, $
               NCOLORS=ncolors, $
               NEGATIVE=negative, $
               TYPE=scaletype, $
               MINTHRESH=sclmin, $
               MAXTHRESH=sclmax, $
               SIGMA=sigma, $
               NOTIFY_OBJ={object:self, method:'XStretch_Notification'}
      
            END
                  
       ELSE  : BEGIN

         event.id -> GetProperty, Value=val
         self._statusbar -> SetProperty, Text= 'Event in IMGWIN'
         ENDCASE

   ENDCASE

END
;*****************************************************************************************************


PRO ImgWin::GUI, menuBar, $
        XWINSIZE=xwinsize, $
        YWINSIZE=ywinsize, $
        BACKGROUND=background, $
        AXES=axes, $
        FULL_RESOLUTION=full_resolution
        
; The purpose of this method is create all the graphical user interface elements.
; Widgets that cause events are named, and the EventHandler method descriminates
; based on the NAME field of the event structure.
   
   ; Create a Quit button in the menu bar.
   fileMenu = OBJ_NEW ('ButtonWidget', menuBar ,  Value='File', /MENU)
   
   ; Add color change and scaling buttons if image is 2D.
   self.theImage -> GetProperty, N_DIMENSIONS=n_dims, XSIZE=imgXsize, YSIZE=imgYsize
   IF n_dims EQ 2 THEN BEGIN
        button = OBJ_NEW('ButtonWidget', fileMenu, Name='IMAGE_COLORS', Value='Change Image Colors')
        button = OBJ_NEW('ButtonWidget', fileMenu, Name='SCALE_IMAGE', Value='Scale Image')
        separator = 1
   ENDIF ELSE separator = 0
   saveasID = OBJ_NEW ('ButtonWidget', fileMenu,  Value='Save Window As...', /MENU, SEPARATOR=separator)
   button = OBJ_NEW ('ButtonWidget', saveasID,  Name='SAVE_WINDOW', Value='JPEG file', UVALUE='JPEG')
   button = OBJ_NEW ('ButtonWidget', saveasID,  Name='SAVE_WINDOW', Value='TIFF File', UVALUE='TIFF')
   button = OBJ_NEW ('ButtonWidget', saveasID,  Name='SAVE_WINDOW', Value='BMP File', UVALUE='BMP')
   button = OBJ_NEW ('ButtonWidget', saveasID,  Name='SAVE_WINDOW', Value='PNG File', UVALUE='PNG')
   button = OBJ_NEW ('ButtonWidget', saveasID,  Name='SAVE_WINDOW', Value='PostScript File', UVALUE='POSTSCRIPT')
   
   properties = OBJ_NEW ('ButtonWidget', fileMenu,  Value='Properties...', /MENU, SEPARATOR=1)
   button = OBJ_NEW ('ButtonWidget', properties,  Name='IMAGE_PROPERTIES', Value='Image Properties')
   IF Keyword_Set(axes) THEN button = OBJ_NEW ('ButtonWidget', properties,  Name='AXES_PROPERTIES', $
        Value='Axes Properties')
 
    exitBttn = OBJ_NEW ('ButtonWidget', fileMenu,  Name='Exit', Value='Exit', /Separator)
   
    IF Keyword_Set(full_resolution) THEN BEGIN
        drawObj = OBJ_NEW ('SelectableDrawWidget', self, XSIZE=imgXsize, YSize=imgYsize, $
             Erase_Window=1, INITIAL_COLOR=background, BUTTON_EVENTS=1, Name='IMGWIN_DRAWWIDGET', $
             /Notify_Realize, MOTION_EVENTS=1, /SCROLL, X_SCROLL_SIZE=xwinsize, Y_SCROLL_SIZE=ywinsize)
    ENDIF ELSE BEGIN
        drawObj = OBJ_NEW ('SelectableDrawWidget', self, XSIZE=xwinsize, YSize=ywinsize, $
             Erase_Window=1, INITIAL_COLOR=background, BUTTON_EVENTS=1, Name='IMGWIN_DRAWWIDGET', $
             /Notify_Realize, MOTION_EVENTS=1)
    ENDELSE
    drawObj -> Add, self.theImage
    self.theDrawWidget = drawObj
   
    ; Create the statusbar.
    self -> CreateStatusBar

    ; Display the entire application in the window.
    self -> Draw, /Center

END
;*****************************************************************************************************


PRO ImgWin::XStretch_Notification, info

; When stretch parameters are changed in XSTRETCH, those parameters are bundled up
; in an info structure, which is passed to this method. Use the info parameters to 
; configure the image.

    @cat_pro_error_handler
    
    ; Pass the new XSTRETCH parameters on to the image.
    IF StrUpCase(info.event_handler) EQ 'XSTRETCH_STRETCHTYPE' THEN BEGIN
        self.theImage -> SetProperty, SCALETYPE=info.type, SCLMIN=info.minThresh, SCLMAX=info.maxThresh, $
               GAMMA=info.gamma, BETA=info.beta, MEAN=info.mean, EXPONENT=info.exponent
    ENDIF ELSE BEGIN
        self.theImage -> SetProperty, SCLMIN=info.minThresh, SCLMAX=info.maxThresh
    ENDELSE
          
    ; Redraw the image in the window
    self.theDrawWidget -> SetWindow
    self.theImage -> Draw
    
END
;*****************************************************************************************************


PRO ImgWin::CLEANUP
    
    @cat_pro_error_handler
    
    Obj_Destroy, self.theImage
    Obj_Destroy, self.theDrawWidget
    Obj_Destroy, self.theAxes
    Obj_Destroy, self._statusbar
    
    self -> TOPLEVELBASE::Cleanup
    
    self -> Report, /Completed
    
END
;*****************************************************************************************************


FUNCTION ImgWin::INIT, image, $
    AXES=axes, $
    BACKGROUND=background, $
    BETA=beta, $
    BOTTOM=bottom, $
    BREWER=brewer, $
    COLOR=color, $
    CTINDEX=ctindex, $
    EXPONENT=exponent, $
    GAMMA=gamma, $
    KEEP_ASPECT=keep_aspect, $
    MEAN=mean, $
    MISSING_COLOR=missing_color, $
    NCOLORS=ncolors, $
    NOINTERPOLATION=nointerp, $
    POSITION=position, $
    SCALETYPE=scaletype, $
    SCLMIN=sclmin, $
    SCLMAX=sclmax, $
    SIGMA=sigma, $
    XRANGE=xrange, $
    XSIZE=xwinsize, $
    XTICKFORMAT=xtickformat, $
    XTITLE=xtitle, $
    YRANGE=yrange, $
    YSIZE=ywinsize, $
    YTICKFORMAT=ytickformat, $
    YTITLE=ytitle, $
    OUTIMAGE=theImage, $   An output keyword.
    _Ref_Extra=extra
    
    ; Catch error handling.
    Catch, theError
    IF theError NE 0 THEN BEGIN
        Catch, /CANCEL
        void = Error_Message()
        RETURN, 0
    ENDIF
    
    ; Must have an image to display. This could be the name of an image file,
    ; or it could be the image itself. Give the user a change to pick an image
    ; file is an image is not specified.
    IF N_Elements(image) EQ 0 THEN BEGIN
       imageFile = Dialog_Pickfile(TITLE='Select image file...')
       IF imageFile EQ "" THEN RETURN, 0
       image = Read_Image(imageFile, r, g, b)
       IF N_Elements(image) EQ 1 THEN Message, 'Image file could not be read with READ_IMAGE.'
       IF N_Elements(r) NE 0 THEN colorPalette = [[r], [g], [b]]
    ENDIF
    
    ; Is this a filename?
    IF Size(image, /TNAME) EQ 'STRING' THEN BEGIN
       imagefile = image
       image = Read_Image(imageFile, r, g, b)
       IF N_Elements(image) EQ 1 THEN Message, 'Image file could not be read with READ_IMAGE.'
       IF N_Elements(r) NE 0 THEN colorPalette = [[r], [g], [b]]
    ENDIF
    
    ; If the image is not a 2D array or true-color image, then return.
    ndims = Size(image, /N_DIMENSIONS)
    IF ndims EQ 3 THEN BEGIN
       index = Where(Size(image, /DIMENSIONS) EQ 3, count)
       IF count LT 0 THEN Message, 'Image does not appear to be a true-color image.'
    ENDIF ELSE BEGIN
       IF ndims LT 2 OR ndims GT 3 THEN Message, 'Image does not appear to be a 2D array.'
    ENDELSE
    
    ; This INIT method simply instantiates a top-level base object with a status bar.
    ok = self->TOPLEVELBASE::INIT(_Extra=extra)
    IF ~ok THEN RETURN, 0
    
    ; Check keywords.
    IF N_Elements(background) EQ 0 THEN background = 'ivory'
    brewer = Keyword_Set(brewer)
    IF N_Elements(color) EQ 0 THEN color = 'charcoal'
    IF N_Elements(ctindex) THEN ctindex = 0
    IF N_Elements(keep_aspect) EQ 0 THEN keep_aspect = 1
    keep_aspect = Keyword_Set(keep_aspect)
    IF N_Elements(nointerp) EQ 0 THEN nointerp = 1
    nointerp = Keyword_Set(nointerp)
    IF N_Elements(position) EQ 0 THEN BEGIN
        IF Keyword_Set(axes) THEN BEGIN
           position = [0.10, 0.075, 0.925, 0.925]
           IF N_Elements(xtitle) NE 0 THEN position[0] = 0.150
           IF N_Elements(ytitle) NE 0 THEN position[1] = 0.125
        ENDIF ELSE BEGIN
           position = [0,0,1,1]
        ENDELSE
    ENDIF
    IF N_Elements(scaletype) EQ 0 $
        THEN scaletype = 0 $
        ELSE BEGIN
             IF Size(scaletype, /TNAME) EQ 'STRING' THEN BEGIN
                  possibleTypes = ['LINEAR', 'GAMMA', 'LOG', 'ASINH', $
                                   'LINEAR 2%', 'SQUARE ROOT', 'EQUALIZATION', 'GAUSSIAN']
                  index = Where(possibleTypes EQ StrUpCase(scaletype), count)
                  IF count EQ 0 THEN Message, 'Unknown scaling type encountered.'
                  scaletype = index
            ENDIF
        ENDELSE

    ; Information about the image.
    dims = Image_Dimensions(image, XSize=xsize, YSize=ysize, TrueIndex=trueindex, $  
       XIndex=xindex, YIndex=yindex)
    imgAspect = Float(ysize) / xsize
    
    ; Check ranges.
    IF N_Elements(xrange) EQ 0 THEN xrange = [0, xsize]
    IF N_Elements(yrange) EQ 0 THEN yrange = [0, ysize]
    
    ; If a colorpalette has not been determined already (from READ_IMAGE), then do it here.
    IF N_Elements(colorpalette) EQ 0 THEN BEGIN
    
        IF N_Elements(ctindex) EQ 0 THEN BEGIN
            ctindex = 0
            TVLCT, rr, gg, bb, /Get
            CTLoad, ctindex, BREWER=brewer
            TVLCT, r, g, b, /GET
            colorPalette = [[r], [g], [b]]
            TVLCT, rr, gg, bb
        ENDIF ELSE BEGIN  
            IF ctindex LT 0 THEN BEGIN
                TVLCT, r, g, b, /GET
                colorPalette = [[r], [g], [b]]
            ENDIF ELSE BEGIN
                TVLCT, rr, gg, bb, /Get
                CTLoad, ctindex, BREWER=brewer
                TVLCT, r, g, b, /GET
                colorPalette = [[r], [g], [b]]
                TVLCT, rr, gg, bb
            ENDELSE   
        ENDELSE
    ENDIF
    
    ; Calculate window size.
    IF (N_Elements(xwinsize) EQ 0) AND (N_Elements(ywinsize) EQ 0) THEN BEGIN
        maxSize = 600 
        maximageSize = Max([xsize, ysize])
        IF 2*maximageSize LT maxsize THEN maxsize = 350 > 2*maximagesize
        IF imgAspect GE 1 THEN BEGIN
            ywinsize = maxsize
            xwinsize = maxsize / imgAspect
        ENDIF ELSE BEGIN
            xwinsize = maxsize
            ywinsize = maxsize * imgAspect
        ENDELSE
    ENDIF ELSE BEGIN
        IF N_Elements(xwinsize) EQ 0 THEN xwinsize = 600
        IF N_Elements(ywinsize) EQ 0 THEN ywinsize = 600
    ENDELSE
    
    ; Create a coordinate system for the image.
    coords = Obj_New('CATCOORD', Name='IMG WIN COORDS OBJECT', XRANGE=xrange, YRANGE=yrange)
    
    ; If AXES are required.
    IF Keyword_Set(axes) THEN BEGIN
       theAxes = Obj_New('IMGAXES', $
            NAME='IMGWIN_AXES', $
            COLOR=color, $
            COORD_OBJECT=coords, $
            POSITION=position, $
            XRANGE=xrange, $
            XTICKFORMAT=xtickformat, $
            XTITLE=xtitle, $
            YRANGE=yrange, $
            YTICKFORMAT=ytickformat, $
            YTITLE=ytitle)
        self.theAxes = theAxes
    ENDIF

    ; Create the image.
    theImage = Obj_New('ScaleImage', image, $
        AXES=theAxes, $
        NAME='IMGWIN_IMAGE', $
        BETA=beta, $
        BOTTOM=bottom, $
        COORD_OBJECT=coords, $
        EXPONENT=exponent, $
        GAMMA=gamma, $
        KEEP_ASPECT=keep_aspect, $
        MEAN=mean, $
        MISSING_COLOR=missing_color, $
        MISSING_VALUE=missing_value, $
        NCOLORS=ncolors, $
        NOINTERPOLATE=nointerp, $
        POSITION=position, $
        SCALETYPE=scaletype, $
        SCLMIN=sclmin, $
        SCLMAX=sclmax, $
        SELECTABLE=1, $
        SIGMA=sigma)
        
    ; Update the colors.
    IF Obj_Valid(theAxes) THEN theImage -> Add, theAxes
    theImage -> GetProperty, COLOR_OBJECT=colors
    colors -> SetProperty, COLORPALETTE=colorPalette, BREWER=Keyword_Set(brewer)
        
    
    ; Store the image.
    self.theImage = theImage
    self.full_resolution = Keyword_Set(full_resolution)
    
    ; Success!
    RETURN, 1
END
;*****************************************************************************************************


PRO ImgWin__Define, class

; The IMGWIN class definition. It is a top-level base object widget.

   compile_opt idl2

   class = { IMGWIN, $
             INHERITS TopLevelBase, $     ; This is an application window.
             theImage: Obj_New(), $       ; The image object to display (a SCALEIMAGE object).
             theDrawWidget: Obj_New(), $  ; The application draw widget.
             theAxes: Obj_New(), $        ; The image axes object (a IMGAXIS object), if required.
             full_resolution: 0B, $       ; Flag that indicates the image is at full-resolution.
             _statusbar: Obj_New()}       ; A status bar for program updates.
END
;*****************************************************************************************************

; This is the driver program for the IMGWIN object.
PRO ImgWin, image, $
    AXES=axes, $
    BACKGROUND=background, $
    BETA=beta, $
    BOTTOM=bottom, $
    BREWER=brewer, $
    COLOR=color, $
    CTINDEX=ctindex, $
    EXPONENT=exponent, $
    FULL_RESOLUTION=full_resolution, $
    GAMMA=gamma, $
    KEEP_ASPECT=keep_aspect, $
    MEAN=mean, $
    MISSING_COLOR=missing_color, $
    NCOLORS=ncolors, $
    NOINTERPOLATION=nointerp, $
    POSITION=position, $
    SCALETYPE=scaletype, $
    SCLMIN=sclmin, $
    SCLMAX=sclmax, $
    SIGMA=sigma, $
    XRANGE=xrange, $
    XSIZE=xwinsize, $
    XTICKFORMAT=xtickformat, $
    XTITLE=xtitle, $
    YRANGE=yrange, $
    YSIZE=ywinsize, $
    YTICKFORMAT=ytickformat, $
    YTITLE=ytitle, $
    OUTIMAGE=outimage  ; An output keyword

   @cat_pro_error_handler
   
   scroll = 0
   IF Keyword_Set(full_resolution) THEN BEGIN
      dims = Image_Dimensions(image, XSIZE=xsize, YSIZE=ysize)
      screenDims = Get_Screen_Size()
      IF xsize GT (screenDims[0] - 8) THEN scroll = 1
      IF ysize GT (screenDims[1] - 95) THEN scroll = 1
   ENDIF

   ; Create the widgets that make up the application. Run it.
   tlb = OBJ_NEW('IMGWIN', image, Column=1, NAME='IMGWIN_TLB', $
        SIZE_EVENTS=1, $
        MBar=menubar, Title='Catalyst Image Window', $
        AXES=axes, $
        BACKGROUND=background, $
        BETA=beta, $
        BOTTOM=bottom, $
        BREWER=brewer, $
        COLOR=color, $
        CTINDEX=ctindex, $
        EXPONENT=exponent, $
        GAMMA=gamma, $
        KEEP_ASPECT=keep_aspect, $
        MEAN=mean, $
        MISSING_COLOR=missing_color, $
        NCOLORS=ncolors, $
        NOINTERPOLATION=nointerp, $
        POSITION=position, $
        SCALETYPE=scaletype, $
        SCLMIN=sclmin, $
        SCLMAX=sclmax, $
        SCROLL=scroll, $
        SIGMA=sigma, $
        XRANGE=xrange, $
        XSIZE=xwinsize, $
        XTICKFORMAT=xtickformat, $
        XTITLE=xtitle, $
        YRANGE=yrange, $
        YSIZE=ywinsize, $
        YTICKFORMAT=ytickformat, $
        YTITLE=ytitle, $
        ; Output keywords
        OUTIMAGE=outimage, $
        _Ref_Extra=extra)

   ; Pass information you will need in the GUI method. I only need it once, so I 
   ; am not storing it in the object.
   IF Obj_Valid(tlb) THEN tlb -> GUI, menubar, $
        XWINSIZE=xwinsize, $
        YWINSIZE=ywinsize, $
        BACKGROUND=background, $
        AXES=axes, $
        FULL_RESOLUTION=full_resolution
        
   
END
;*****************************************************************************************************
