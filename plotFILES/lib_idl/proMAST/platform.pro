;******************************************************************************
;+
;*NAME:
;
;       PLATFORM
;
;*CLASS:
;
;       system dependancies
;
;*CATEGORY:
;
;*PURPOSE:
;
;       Set system dependant items such as system commands that will be
;	spawned, internal data formats, etc.
;
;*CALLING SEQUENCE:
;
;       PLATFORM,DUMMY,ncopy=ncopy,icopy=icopy,copymes=copymes,   $
;          ndelete=ndelete,idelete=idelete,pdelete=pdelete,npurge=npurge,   $
;          rename=rename,page=page,cat=cat,submit=submit,awk=awk,node=node,   $
;          domain=domain,idfmt=idfmt,rms=rms,keyst=keyst,sysos=sysos,   $
;          mail=mail,printcmd=printcmd,mainprt=mainprt,secprt=secprt,   $
;          cprint=cprint,dirlis=dirlis,find=find,security=security,   $
;          quota=quota,syntax=syntax,font=font,syslogin=syslogin,   $
;          trmldev=trmldev,nulldev=nulldev,tapectrl=tapectrl,   $
;          tapecopy=tapecopy,tapetar=tapetar,tapedrive=tapedrive, $
;          tmploc=tmploc
;
;*PARAMETERS:
;
;       DUMMY   (REQ) (I) () ()
;               A dummy variable that will allow the calling sequence to be
;               printed.
;
;       NCOPY	(KEY) (O) (0) (S)
;               The copy command without conformation (non-interactive).
;
;       ICOPY	(KEY) (O) (0) (S)
;               The copy command with conformation (interactive).
;
;       COPYMES (KEY) (O) (0) (S)
;               Message related to copy command.
;                  copymes.present - present tense of the word.
;                  copymes.past - past tense of the word.
;
;       NDELETE	(KEY) (O) (0) (S)
;               The delete command without conformation (non-interactive).
;
;       IDELETE (KEY) (O) (0) (S)
;               The delete command with conformation (interactive).
;
;	PDELETE	(KEY) (O) (0) (S)
;		Print file and then delete file.
;
;       NPURGE  (KEY) (O) (0) (S)
;               The purge command without conformation (non-interactive).
;       
;	RENAME  (KEY) (O) (0) (S)
;		To give a file a different name (rename or move).
;
;       PAGE	(KEY) (O) (0) (S)
;               The page command for displaying a file.
;
;       CAT     (KEY) (O) (0) (S)
;               The cat or type command.
;
;       SUBMIT  (KEY) (O) (0) (S)
;               The command to submit a batch job.
;                  submit.cmd - the submit command
;                  submit.notify - the notify option
;                  submit.noprint - the noprint option
;
;       AWK     (KEY) (O) (0) (S)
;               The awk command information.
;                  awk.cmd - the command
;                  awk.file - the file of awk commands
;
;       NODE    (KEY) (O) (0) (S)
;               The system's name.  On VMS systems, obtained via
;               getenv('NODE').  On UNIX systems, obtained via getenv('HOST').
;               If none obtained, then set to ''.
;
;       DOMAIN  (KEY) (O) (0) (S)
;               The Internet domain - when applicable.  Needs to be changed
;               by the person responsible for the software if the correct
;               domain is NOT '.gsfc.nasa.gov'.
;
;       IDFMT	(KEY) (I/O) (0) (S)
;               The internal data format.  Types are:
;                  'ieeebig' - IEEE big-endian format (SunOS/sparc, FITS)
;                  'ieeelittle' - IEEE little-endian format (Ultrix/mipsel,
;                                 DOS/386, VMS/mipsel)
;                  'vaxformat' - byte, signed word and longword integers,
;                                f-floating, d-floating (Ultrix/VAX, VMS/VAX)
;		The variable must be sent to ' ' before set to PLATFORM for
;		this parameter to be returned.
;
;       RMS     (KEY) (O) (0) (I)
;               Record Management Services.  For use with some (CON_RDAF and
;               CONVERT) of the internal data format conversion procedures.
;               Set to 1 when the system is VMS or Ultrix and the architecture
;               is VAX or Mipsel.  For all other cases, set to 0.  Used to
;               consider Ultrix/VAX to VMS/VAX and Ultrix/Mipsel and VMS/Mipsel
;               "translations".
;
;	KEYST	(KEY) (O) (0) (B)
;		The keystroke value for carriage return.
;
;	SYSOS	(KEY) (O) (0) (S)
;		The value of !version.os.  Main uses include
;		- calling sequences on a non-VMS machine will that VMS only
;		  keywords will cause to crash
;		- record lengths
;               - case sensitivity
;               - one more entry in listing file (non-VMS systems)
;               - !p.font setting
;               - state the operating system
;               - most of IUECOPY
;
;       MAIL	(KEY) (O) (0) (S)
;               The mail command.
;                  mail.command - mail command itself with the subject
;                                 parameter.
;                  mail.noself - option to add to surpress cc to the current
;                                account.
;                  mail.redirect - allows a file to be sent to the address
;                                  (Unix).
;                  mail.notallow - characters not allowed in subject part of
;                                  mail command
;                  mail.ndads - email address for the request account at NSSDC
;                              (for IUEFX).
;		   mail.usrreq - email account at Greentec (for USRREQ).
;
;	PRINTCMD (KEY) (O) (0) (S)
;		The default print command.
;
;       MAINPRT (KEY) (O) (0) (S)
;               Main printer information.
;                  mainprt.type - the type (such as 'ps')
;                  mainprt.dev - device commands
;                  mainprt.com - the command
;
;       SECPRT  (KEY) (O) (0) (S)
;               If there is a second printer, its information.
;                  secprt.type - the type (such as 'ps')
;                  secprt.dev - device commands
;                  secprt.com - the command
;
;	CPRINT	(KEY) (O) (0) (S)
;		The color print command.
;
;	DIRLIS	(KEY) (O) (0) (S)
;		Directory listing.
;                  dirlis.cmd - the basic directoty listing command.
;                  dirlis.sec - the option to see security settings
;                  dirlis.col - listing in one column
;                  dirlis.output - output file
;                  dirlis.nohead - only filenames (no extra "stuff" printed)
;                  dirlis.since - limit output to ust files changed in past
;                                 24 hours
;
;       FIND    (KEY) (O) (0) (S)
;               The find command.
;                   find.cmd - the command
;                   find.ctime - how many days since file changed
;                   find.prune - prune search tree
;                   find.type - type of file
;                   find.name - filename to be searched for
;                   find.print - output filenames found
;
;	SECURITY (KEY) (O) (0) (S)
;		Reset secuity permissions on a file.
;                  security.set - decline other users from readi write, or
;                                 execute privilege.
;		   ecurity.unset - allow everyone read permission and group
;		                      execute privilege.
;
;	QUOTA	(KEY) (O) (0) (S)
;		Show disk quota command.
;
;       SYNTAX  (KEY) (O) (0) (S)
;               Structure of different syntax notations.
;                  syntax.redirect - redirect output to a file
;                  syntax.addon - attached information to end
;                  syntax.pipe - pipe output to another command
;                  syntax.background - put process in background
;                  syntax.execute - execute character
;                  syntax.allver - all versions (when versions available)
;                  syntax.disksep - separates disk name from path
;                  syntax.startpath - separates the disk from the start of the
;                                     path
;                  syntax.endpath - separates end of path from filename
;                  syntax.versep - separates filename from version number
;                  syntax.extsep - separates filename extension from the rest
;                                  of the file name
;                  syntax.listsep - separate names in a list of filenames
;                  syntax.pdelim - separate list of files to print command
;                  syntax.pathlist - separate directory names in !path
;
;       FONT    (KEY) (O) (0) (S)
;               The font available for fixed width.
;
;       SYSLOGIN  (KEY) (O) (0) (S)
;               The system login logical.
;
;	TRMLDEV	(KEY) (O) (0) (S)
;		Terminal device name.
;
;       NULLDEV (KEY) (O) (0) (S)
;               The null device (garbage maybe piped to it).
;
;       TAPECTRL (KEY) (O) (0) (S)
;               Magnetic tape control for skipping files - Unix and Ultrix
;               systems.  tape.device = '' is assumed to be a VMS or VMS like
;               system.
;                  tapectrl.device - the device directory.
;                  tapectrl.mt - the tape control command.
;                  tapectrl.back - back up on the tape.
;                  tapectrl.forward - go forward on the tape.
;                  tapectrl.bandf - back and then forward one.
;                  tapectrl.fspace - forward space on a tape.
;                  tapectrl.bspace - back space on a tape.
;                  tapectrl.status - print status information about the tape
;                                    unit.
;                  tapectrl.rewind - rewind the tape.
;                  tapectrl.weof - write EOF marker
;                  tapectrl.eom - goto end of recorded media
;                  tapectrl.eom - goto end of recorded media
;
;	TAPECOPY (KEY) (O) (0) (S)
;               Tape copy commands.
;                  tapecopy.cmd - command.
;                  tapecopy.infile - input file marker.
;                  tapecopy.outfile - output file marker.
;                  tapecopy.ibs - input block size marker.
;                  tapecopy.obs - output size marker.
;                  tapecopy.num - number of files to copy before ending marker.
;                  tapecopy.qic - block size for QIC tapes (must be 512).
;                  tapecopy.nonqic - block size for non-QIC tapes.
;                  tapecopy.conv - conversion to take place
;                  tapecopy.sync - with tapecopy.conv, pad every input record
;                                  to ibs
;
;       TAPETAR    (KEY) (O) (0) (I)
;                  Equals 1 if tar is available.  Equals 0 is tar is not
;                  available.
;
;       TAPEDRIVE  (KEY) (O) (1) (S)
;                  The array of tape drive(s) available.  If no tape drive are
;                  available, tapedrive(0) = ''.  The devices listed are
;                  specific to the LASP/GSFC VAX Cluster that the IUE VAX is a
;                  member of and IUE SUN 4/280, iuesn1.  Other facilities
;                  should modify this value to describe their own available
;                  tape drive(s).
;
;       TMPLOC     (KEY) (O) (1) (S)
;                  input string specifying a path for temporary scratch files.
;
;*EXAMPLES:
;
;*SYSTEM VARIABLES USED:
;
;	!version.os
;	!d.name
;
;*INTERACTIVE INPUT:
;
;	none
;
;*SUBROUTINES CALLED:
;
;       PARCHECK
;       INTDAFMT
;       VMSPLAT
;       UIXPLAT
;       MACPLAT
;       DOSPLAT
;
;*FILES USED:
;
;       none
;
;*SIDE EFFECTS:
;
;*RESTRICTIONS:
;
;	In INTDAFMT, if !version.arch is not one of the given values, the
;       procedure tries to determine the internal data format.
;
;	If !version.os is not one of the given values, the procedure defaults
;	to the Unix values.
;
;*NOTES:
;
;       Please see the file platform.txt in the iuerdaf/manual subdirectory
;       for information on how to customizing this procedure.
;
;	The VMS branch assumes Decnet access for mail.  Examples for Internet
;	access are given.  The other branches assume Internet access.  Please
;       consult your system manager if you are unsure.
;
;	Additional systems may be added by contacting IUEDAC staff.  However,
;       please note that they will be UNsupported.  Please include the
;       information requested in the platform.txt file.
;
;*PROCEDURE:
;
;       The HOST and NODE environment variables are checked to determine the
;       system's name.
;
;	If the IDFMT keyword is set (why it must be set to ' ' before being
;	requested), the internal data format is determined in INTDAFMT.  If the
;	!version.arch is not one of the know values, then a test is used to
;	determine which format it is.  If it is not IEEE big-endian, IEEE
;	little-endian, or VAX formats, then the procedure retalls.  If the
;       system is VMS or Ultrix with an architectue of VAX or Mipsel, then
;       the RMS (Record Management Services) keyword is set to 1.  For all
;       other cases, it equals 0.
;
;	The other keywords are set based on !version.os.  If !version.os is not
;       one of the know values, then the Unix version is used as the default.
;
;*I_HELP  nn:
;
;*MODIFICATION HISTORY:
;
;	 6 Jun 94  PJL  wrote
;        9 Jun 94  PJL  added integer check and RMS keyword to internal data
;                       format section
;	10 Jun 94  PJL  changed mtskip to tape and added tags; added tapecopy
;       13 Jun 94  PJL  add tags to tapecopy;  add tapedrive, node, syntax, and
;                       nulldev
;       14 Jun 94  PJL  add tags to syntax;  added submit, find, syslogin, awk,
;                       mainprt, secprt
;       16 Jun 94  PJL  add npurge, ibmallow, and tapetar
;       24 Jun 94  PJL  added domain
;        5 Jul 94  PJL  rearrange some and add addition documentation
;       11 Jul 94  PJL  correct error in integer section of internal data
;                       format
;       25 Jul 94  PJL  correct typo
;       26 Aug 94  RWT  change "tt:/fort" to "sys$output" for vms systems
;	22 Sep 94  PJL  broke off non-installion site dependant items into
;                       VMSPLAT, UOXPLAT, DOSPLAT, MACPLAT, and INTDAFMT
;       26 Sep 94  PJL  removed device 5 - MSA0  (it has been physically
;                       removed); combined ultrix, sunos, and else cases into
;                       else case; added checks on node name for iuesn1
;                       specifics
;        6 Oct 94  PJL  removed ibmallow keyword
;        7 Oct 94  PJL  mail.command, cprint set to 'NA' when not available
;       27 Oct 94  PJL  added SPUD 8mm tape drives
;       28 Nov 94  PJL  updated VMS tape list for LASP cluster; added font
;       20 Mar 95  RWT  add tmploc keyword
;       12 May 95  PJL  changed else mail command to /usr/ucb/mail to
;			allow the '-s' command
;       19 May 95  RWT  set mainprt.type to 'PS' for Macs and PCs
;       30 May 95  RWT  set mainprt.device for Macs and Pcs
;       24 Jul 95  RWT  modify for IDL 4.0 !version changes for PCs
;        7 May 96  RWT  point to NDADS rather than NDADSA
;	 2 Dec 96  PJL  modified for OSF system
;       18 jun 97  jrc  branch in unix mail command for osf system
;-
;******************************************************************************
 pro platform,dummy,ncopy=ncopy,icopy=icopy,copymes=copymes,ndelete=ndelete, $
        idelete=idelete,pdelete=pdelete,npurge=npurge,rename=rename,page=page, $
        cat=cat,submit=submit,awk=awk,node=node,domain=domain,idfmt=idfmt,   $
        rms=rms,keyst=keyst,sysos=sysos,mail=mail,printcmd=printcmd,   $
        mainprt=mainprt,secprt=secprt,cprint=cprint,dirlis=dirlis,find=find,   $
        security=security,quota=quota,syntax=syntax,font=font,   $
        syslogin=syslogin,trmldev=trmldev,nulldev=nulldev,tapectrl=tapectrl,  $
        tapecopy=tapecopy,tapetar=tapetar,tapedrive=tapedrive,tmploc=tmploc
;
 npar = n_params(0)
 if (npar eq 0) then begin
    print,'PLATFORM,DUMMY,ncopy=ncopy,icopy=icopy,copymes=copymes,' +    $
       'ndelete=ndelete,   $'
    print,'   npurge=npurge,idelete=idelete,pdelete=pdelete,rename=rename,' + $
       'page=page,   $'
    print,'   cat=cat,submit=submit,awk=awk,node=node,domain=domain,' +   $
       'idfmt=idfmt,rms=rms,  $'
    print,'   keyst=keyst,sysos=sysos,mail=mail,printcmd=printcmd,' +   $
       'mainprt=mainprt,   $'
    print,'   secprt=secprt,cprint=cprint,dirlis=dirlis,find=find,' +   $
       'security=security,   $'
    print,'   quota=quota,syntax=syntax,font=font,syslogin=syslogin,' +   $
       'trmldev=trmldev,   $'
    print,'   nulldev=nulldev,tapectrl=tapectrl,tapecopy=tapecopy,' +   $
       'tapetar=tapetar,   $'
    print,'   tapedrive=tapedrive,tmploc=tmploc'
    retall
 endif  ; npar eq 0
 parcheck,npar,1,'PLATFORM'
;
;  determine the node name
;
 temp = strlowcase(strtrim(getenv('HOST'),2))
 if (temp ne '') then node = temp else begin
    temp = strlowcase(strtrim(getenv('NODE'),2))
    if (temp ne '') then node = temp else node = ''
 endelse  ; temp ne ''
;
;  operating system
;
 sysos = strlowcase(!version.os)
 if (sysos eq 'win32') then sysos = 'windows'
;
;  internal data format
;
 if (keyword_set(idfmt)) then intdafmt,idfmt,rms
;
; check print queue logical - for a qms printer
;
 pque = getenv('LASER')     ; check for qms print que
;
;  commands determined by the operating system - see the prolog for information
;
 case sysos of
;
;  VMS types
;    { vax vms 3.5.1}
;    { alpha vms 3.5.1}
;
    'vms':  begin
;
;  mail - assuming DECnet
;
               mail = {mail,command:'mail/subject=',noself:'/noself ',   $
                       redirect:' ',notallow:'',ndads:' NDADS::ARCHIVES ',   $
                       usrreq:' IUEGTC::USRREQ '}
;
;  For Internet replace the ndads and usrreq entries above with something like:
;                       ndads:' SMTP%"""archives@ndads.gsfc.nasa.gov""" '
;                       usrreq:' SMTP%"""usrreq@iuegtc.gsfc.nasa.gov""" '
;  for systems using MultiNet or
;                       ndads:' WINS%"""archives@ndads.gsfc.nasa.gov""" '
;                       usrreq:' WINS%"""usrreq@iuegtc.gsfc.nasa.gov""" '
;  for systems using Wollongong.  Consult your system manager for assistance.
;
;  domain not needed for DECnet - domain needed for Internet
;
               domain = ''
;
;  printers
;     mainprt is the main laser printer - use IUELIB to set or reset
;     secprt is the second laser printer - use QMS to change to
;
               psque = getenv('PSLASER')  ; check for postscript print que
               if (psque ne '') then begin
;
;  postscript printer (default)
;
                  mainprt = {mainprt,type:'ps',dev:'/landscape,filename=',   $
                             com:'print/que=PSLASER'}
                  if (pque ne '') then begin
                     if (strlowcase(!d.name) eq 'tek') then secprt =   $
                        {secprt,type:'tek',dev:'/tty,filename=',   $
                         com:'print/que=LASER'}
                     if (strlowcase(!d.name) ne 'tek') then secprt =    $
                        {secprt,type:'tek',dev:'filename=',   $
                         com:'print/que=LASER'}
                  endif else secprt = {secprt,type:'',dev:'',com:''}
               endif else begin
;
;  QMS printer with tek or non-tek screen:
;
                  if (pque ne '') then begin
                     if (strlowcase(!d.name) eq 'tek') then mainprt =   $
                        {mainprt,type:'tek',dev:'/tty,filename=',   $
                         com:'print/que=LASER'}
                     if (strlowcase(!d.name) ne 'tek') then mainprt =   $
                        {mainprt,type:'tek',dev:'filename=',   $
                         com:'print/que=LASER'}
                  endif else begin
                     mainprt = {mainprt,type:'',dev:'',com:''}
                     secprt = {secprt,type:'',dev:'',com:''}
                  endelse  ; pque ne ''
               endelse  ; psque ne ''
;
;  color printer command
;
               cprint = 'print/que=iue$teksd '
;
;  IUE tape drives
;
               tapedrive = ['unit     name',   $
                            ' 1     STARS$MUA0',   $
                            ' 2     STARS$MUB0',   $
                            ' 3     SPUD$MKB0  (8mm - only accessible ' +   $
                                                'from SPUD)',  $
                            ' 4     $2$MUA0',   $
                            ' 5     SPUD$MKB100  (8mm - only accessible ' +   $
                                                  'from SPUD)', $
                            ' 6     MSB0',   $
                            ' 7     $2$MUA1',   $
                            ' 8     $2$MUA2']
;
;  CASA VAX tape drives
;
;;;               tapedrive = ['   From system HYADES - ',   $
;;;               '     unit     name ',   $
;;;               '      1       MUB1 - 8mm  Exabyte drive, 2.3 GB ',   $
;;;               '   From system TAURUS - ',   $
;;;               '     unit     name ',   $
;;;               '      0       MSA0 - 1/2" tape drive, 800 to 1600 bpi ',   $
;;;               '      1       MSB0 - 1/2" tape drive, 800 to 1600 bpi ',   $
;;;               '   From system CYGNUS - ',   $
;;;               '     unit     name ',   $
;;;               '      0       MSA0 - 1/2" tape drive, 800 to 6250 bpi ',   $
;;;               '      1       MSB0 - 1/2" tape drive, 800 to 6250 bpi ']
;
;  the following should be correct for all VMS systems
;
;  tar command is usually not available on VMS systems
;
               tapetar = 0
;
;  other vms commands
;
               vmsplat,ncopy,icopy,copymes,ndelete,idelete,pdelete,   $
                  npurge,rename,page,cat,submit,awk,keyst,printcmd,dirlis,   $
                  find,security,quota,syntax,font,syslogin,trmldev,nulldev,   $
                  tapectrl,tapecopy,tmploc
            end  ; !version.os eq 'vms'
;
;  DOS/Windows PC types
;    { 3.1 windows 3.5.1}
;    { x86 win32   4.0   windows}
;
    'windows':  begin
;
;  assuming Internet
;
                   mail = {mail,command:'NA',noself:'',redirect:'',   $
                           notallow:'',   $
                           ndads:' archives@ndads.gsfc.nasa.gov ',   $
                           usrreq:' usrreq@iuegtc.gsfc.nasa.gov '}
;
;  for Internet, domain needed
;
                   domain = '.gsfc.nasa.gov'
;
;  printers
;     mainprt is the main laser printer - use IUELIB to set or reset
;     secprt is the second laser printer - use QMS to change to
;
                   mainprt = {mainprt,type:'ps',dev:'/landscape,filename=', $
                              com:''}
                   secprt = {secprt,type:'',dev:'',com:''}
;
;  color printer command
;
                   cprint = 'NA'
;
;  assuming no tape drives available - add as in sections above for given
;  system
;
                   tapedrive = ['']
;
;  the following should be correct for all DOS/Windows systems
;
;  tar command is usually not available
;
                   tapetar = 0
;
;
;  other dos/windows commands
;
                   dosplat,ncopy,icopy,copymes,ndelete,idelete,pdelete,   $
                      npurge,rename,page,cat,submit,awk,keyst,printcmd,   $
                      dirlis,find,security,quota,syntax,font,syslogin,   $
                      trmldev,nulldev,tapectrl,tapecopy,tmploc
                end  ; !version.os eq ('windows' or 'win32')
;
;  Macs
;    { PowerMac MacOS 3.6.1b
;    { Macintosh MacOS 3.5.1}
;
   'macos':  begin
;
;  assuming Internet
;
                mail = {mail,command:'NA',noself:'',redirect:'',notallow:'',  $
                        ndads:' archives@ndads.gsfc.nasa.gov ',   $
                        usrreq:' usrreq@iuegtc.gsfc.nasa.gov '}
;
;  for Internet, domain needed
;
                domain = '.gsfc.nasa.gov'
;
;  printers
;     mainprt is the main laser printer - use IUELIB to set or reset
;     secprt is the second laser printer - use QMS to change to
;
                mainprt = {mainprt,type:'ps',dev:'/landscape,filename=',com:''}
                secprt = {secprt,type:'',dev:'',com:''}
;
;  color printer command
;
                cprint = 'NA'
;
;  assuming no tape drives available - add as in sections above for given
;  system
;
                tapedrive = ['']
;
;  the following should be correct for all Mac systems
;
;  tar command is usually available
;
                tapetar = 0
;
;  if applescipt is available, applescript should be set to 1, else to 0
;
                applescript = 0
;
;
;  other Mac commands
;
                macplat,applescript,ncopy,icopy,copymes,ndelete,idelete,   $
                   pdelete,npurge,rename,page,cat,submit,awk,keyst,printcmd,  $
                   dirlis,find,security,quota,syntax,font,syslogin,trmldev,   $
                   nulldev,tapectrl,tapecopy,tmploc
             end  ; !version.os eq 'macos'
;
;  Ultrix types
;    { mipsel ultrix 3.1.0}
;  Unix types
;    { sparc sunos 3.5.1}
;  others - I have no idea - default to unix type until find out otherwise
;    { hp_pa hp-ux 3.5.1}  ;hp unix
;    { mipseb IRIX 3.5.1}  ;sgi
;    { alpha OSF 3.5.1}
;    { ibmr2 AIX 3.5.1}
;
    else:  begin
;
;  assuming Internet
;
              if strlowcase(!version.os) eq 'osf' then $
;
;  different mail command for osf systems
;
                 mail = {mail,command:'/usr/bin/Mail -s ',noself:' ',   $
                         redirect:' < ',notallow:'$',   $
                         ndads:' archives@ndads.gsfc.nasa.gov ',   $
                         usrreq:' usrreq@iuegtc.gsfc.nasa.gov '}   $
              else   $
;
;  use this mail command instead
;
                 mail = {mail,command:'/usr/ucb/mail -s ',noself:' ',   $
                         redirect:' < ',notallow:'$',   $
                         ndads:' archives@ndads.gsfc.nasa.gov ',   $
                         usrreq:' usrreq@iuegtc.gsfc.nasa.gov '}
;
;  for Internet, domain needed
;
              domain = '.gsfc.nasa.gov'
;
;  printers
;     mainprt is the main laser printer - use IUELIB to set or reset
;     secprt is the second laser printer - use QMS to change to
;
              psque = getenv('PRINTER')   ; check for postscript print que
              if (psque ne '') then begin
;
; postscript printer (default):
;
                 mainprt = {mainprt,type:'ps',dev:'/landscape,filename=', $
                            com:'lpr -h'}
;
;  for a second laser printer - assuming QMS printer with tek or non-tek screen
;
                 if (pque ne '') then begin
                     if (strlowcase(!d.name) eq 'tek') then secprt =   $
                        {secprt,type:'tek',dev:'/tty,filename=',   $
                         com:'lpr -h -Pqms'}
                     if (strlowcase(!d.name) ne 'tek') then secprt =   $
                        {secprt,type:'tek',dev:'filename=',   $
                         com:'lpr -h -Pqms'}
;
;  no second laser printer
;
                 endif else secprt = {secprt,type:'',dev:'',com:''}
              endif else begin
;
; QMS printer with tek or non-tek screen:
;
                 if (pque ne '') then begin
                    if (strlowcase(!d.name) eq 'tek') then mainprt =    $
                       {mainprt,type:'tek',dev:'/tty,filename=',   $
                        com:'lpr -h -Pqms'}
                    if (strlowcase(!d.name) ne 'tek') then mainprt =    $
                       {mainprt,type:'tek',dev:'filename=',   $
                        com:'lpr -h -Pqms'}
                 endif else begin
;
;  no primary nor secondary laser printers
;
                    mainprt = {mainprt,type:'',dev:'',com:''}
                    secprt = {secprt,type:'',dev:'',com:''}
                 endelse  ; pque ne ''
              endelse  ; psque ne ''
;
;  for iuesn1 - postscript printer in the trailer and color printer available
;
              if (node eq 'iuesn1') then secprt = {secprt,type:'ps',   $
                 dev:'/landscape,filename=',com:'lpr -h -Ptps'}
;
;  iuesn1 color printer command
;
              if (node eq 'iuesn1') then cprint = 'lpr -h -Pphaser ' else   $
;
;  assume no color printer available
;
                 cprint = 'NA'
;
;  for ndads1 - postscript printers in Buildings 26 and 28 
;     color printer available
;
              if (node eq 'ndads1') then begin
                 mainprt = {mainprt,type:'ps',dev:'/landscape,filename=', $
                            com:'lpr -Plp1 -h'}
                 secprt = {secprt,type:'ps',   $
                    dev:'/landscape,filename=',com:'lpr -h -Ptps'}
              endif  ; node eq 'ndads1'
;
;  iuesn1 tape drives
;
              if (node eq 'iuesn1') then tapedrive = [   $
                 '     nrmt0 - 1/2" tape drive, 1600 bpi, no auto rewind',  $
                 '     nrmt8 - 1/2" tape drive, 6250 bpi, no auto rewind',  $
                 '     nrst0 - 1/4" cartridge tape drive, QUIC 11 format',  $
                 '     nrst8 - 1/4" cartridge tape drive, QUIC 24 format']  $
;
;  assuming no tape drives available - add as in sections above for given
;  system
;
              else tapedrive = ['']
;
;  ndads1 tape drive
;
              if (node eq 'ndads1') then tapedrive = [   $
                 '     nrmt0l - 8mm tape drive, no auto rewind']
;
;  the following should be correct for all Ultrix, OSF, and SunOS systems
;
;  tar command is usually available
;
              tapetar = 1
;
;  other commands
;
              uixplat,ncopy,icopy,copymes,ndelete,idelete,pdelete,npurge,   $
                 rename,page,cat,submit,awk,keyst,printcmd,dirlis,find,   $
                 security,quota,syntax,font,syslogin,trmldev,nulldev,   $
                 tapectrl,tapecopy,tmploc
           end  ; !version.os - else
;
 endcase  ; !version.os
;
 return
 end  ; platform
