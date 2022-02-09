;******************************************************************************
;+
;*NAME:
;
;       UIXPLAT
;
;*CLASS:
;
;       system dependancies - Unix and Ultrix types
;
;*CATEGORY:
;
;*PURPOSE:
;
;       Set system dependant items such as system commands that will be
;	spawned, internal data formats, etc.
;       Unix and Ultrix types - such as { sparc sunos 3.5.1} and
;       { mipsel ultrix 3.1.0}.  ALso the default.
;
;*CALLING SEQUENCE:
;
;       UIXPLAT,NCOPY,ICOPY,COPYMES,NDELETE,IDELETE,PDELETE,NPURGE,RENAME,  $
;	   PAGE,CAT,SUBMIT,AWK,KEYST,PRINTCMD,DIRLIS,FIND,SECURITY,   $
;          QUOTA,SYNTAX,FONT,SYSLOGIN,TRMLDEV,NULLDEV,TAPECTRL,TAPECOPY,TMPLOC
;
;*PARAMETERS:
;
;       NCOPY	(REQ) (O) (0) (S)
;               The copy command without conformation (non-interactive).
;
;       ICOPY	(REQ) (O) (0) (S)
;               The copy command with conformation (interactive).
;
;       COPYMES (REQ) (O) (0) (S)
;               Message related to copy command.
;                  copymes.present - present tense of the word.
;                  copymes.past - past tense of the word.
;
;       NDELETE	(REQ) (O) (0) (S)
;               The delete command without conformation (non-interactive).
;
;       IDELETE (REQ) (O) (0) (S)
;               The delete command with conformation (interactive).
;
;	PDELETE	(REQ) (O) (0) (S)
;		Print file and then delete file.
;
;       NPURGE  (REQ) (O) (0) (S)
;               The purge command without conformation (non-interactive).
;       
;	RENAME  (REQ) (O) (0) (S)
;		To give a file a different name (rename or move).
;
;       PAGE	(REQ) (O) (0) (S)
;               The page command for displaying a file.
;
;       CAT     (REQ) (O) (0) (S)
;               The cat or type command.
;
;       SUBMIT  (REQ) (O) (0) (S)
;               The command to submit a batch job.
;                  submit.cmd - the submit command
;                  submit.notify - the notify option
;                  submit.noprint - the noprint option
;
;       AWK     (REQ) (O) (0) (S)
;               The awk command information.
;                  awk.cmd - the command
;                  awk.file - the file of awk commands
;
;	KEYST	(REQ) (O) (0) (B)
;		The keystroke value for carriage return.
;
;	PRINTCMD (REQ) (O) (0) (S)
;		The default print command.
;
;	DIRLIS	(REQ) (O) (0) (S)
;		Directory listing.
;                  dirlis.cmd - the basic directoty listing command.
;                  dirlis.sec - the option to see security settings
;                  dirlis.col - listing in one column
;                  dirlis.output - output file
;                  dirlis.nohead - only filenames (no extra "stuff" printed)
;                  dirlis.since - limit output to ust files changed in past
;                                 24 hours
;
;       FIND    (REQ) (O) (0) (S)
;               The find command.
;                   find.cmd - the command
;                   find.ctime - how many days since file changed
;                   find.prune - prune search tree
;                   find.type - type of file
;                   find.name - filename to be searched for
;                   find.print - output filenames found
;
;	SECURITY (REQ) (O) (0) (S)
;		Reset secuity permissions on a file.
;                  security.set - decline other users from readi write, or
;                                 execute privilege.
;		   ecurity.unset - allow everyone read permission and group
;		                      execute privilege.
;
;	QUOTA	(REQ) (O) (0) (S)
;		Show disk quota command.
;
;       SYNTAX  (REQ) (O) (0) (S)
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
;                  syntax.midpath - separates the directory and subdirectory
;                  syntax.endpath - separates end of path from filename
;                  syntax.versep - separates filename from version number
;                  syntax.extsep - separates filename extension from the rest
;                                  of the file name
;                  syntax.listsep - separate names in a list of filenames
;                  syntax.pdelim - separate list of files to print command
;                  syntax.pathlist - separate directory names in !path
;
;       FONT    (REQ) (O) (0) (S)
;               The fixed width font available.
;
;       SYSLOGIN  (REQ) (O) (0) (S)
;               The system login logical.
;
;	TRMLDEV	(REQ) (O) (0) (S)
;		Terminal device name.
;
;       NULLDEV (REQ) (O) (0) (S)
;               The null device (garbage maybe piped to it).
;
;       TAPECTRL (REQ) (O) (0) (S)
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
;	TAPECOPY (REQ) (O) (0) (S)
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
;	TMPLOC (REQ) (O) (0) (S)
;                  string describing the path for storing temporary scratch
;                  files (set to '/tmp' for unix)
;
;*EXAMPLES:
;
;*SYSTEM VARIABLES USED:
;
;	!version.os
;
;*INTERACTIVE INPUT:
;
;	none
;
;*SUBROUTINES CALLED:
;
;       PARCHECK
;
;*FILES USED:
;
;       none
;
;*SIDE EFFECTS:
;
;*RESTRICTIONS:
;
;*NOTES:
;
;       Please see the file platform.txt in the iuerdaf/manual subdirectory
;       for information on how to customizing this procedure.
;
;	Additional systems may be added by contacting IUEDAC staff.  However,
;       please note that they will be UNsupported.  Please include the
;       information requested in the platform.txt file.
;
;*PROCEDURE:
;
;	Keywords are set based on !version.os.
;
;*I_HELP  nn:
;
;*MODIFICATION HISTORY:
;
;       22 Sep 94  PJL  broke out from PLATFORM
;       11 Oct 94  PJL  changed npurge, submit, and quota from '' to 'NA'
;       20 Oct 94  PJL  added syntax.midpath
;       28 Nov 94  PJL  added font
;       20 Mar 95  RWT  add tmploc parameter
;       15 Apr 96  jrc  correct error for dirlis.col to force single column
;                       output
;	 2 Dec 96  PJL  added osf section - WARNING - NOT TESTED
;
;-
;******************************************************************************
 pro uixplat,ncopy,icopy,copymes,ndelete,idelete,pdelete,npurge,rename,page, $
        cat,submit,awk,keyst,printcmd,dirlis,find,security,quota,syntax,font, $
        syslogin,trmldev,nulldev,tapectrl,tapecopy,tmploc
;
 npar = n_params(0)
 if (npar eq 0) then begin
    print,'UIXPLAT,NCOPY,ICOPY,COPYMES,NDELETE,IDELETE,PDELETE,NPURGE,' +   $
       'RENAME,  $'
    print,'   PAGE,CAT,SUBMIT,AWK,KEYST,PRINTCMD,DIRLIS,FIND,SECURITY,   $'
    print,'   QUOTA,SYNTAX,FONT,SYSLOGIN,TRMLDEV,NULLDEV,TAPECTRL,TAPECOPY, $'
    print,'   TMPLOC'
    retall
 endif  ; npar eq 0
 parcheck,npar,24,'UIXPLAT'
;
;  commands determined by the operating system - see the prolog for information
;  defaults to these if !version.os not any of the others
;
;  the following should be correct for all Unix and Ultrix type systems
;
 ncopy = 'echo y | cp -i '
 icopy = 'cp -i '
 copymes = {copymes,present:' overwrite ',past:' overwritten '}
 ndelete = 'echo y | rm -i '
 idelete = 'rm -i '
 pdelete = 'lpr -r '
 npurge = 'NA'
 rename = 'mv '
 page = 'page '
 cat = 'cat '
 submit = {submit,cmd:'',notify:'',noprint:'',delete:''}
 awk = {awk,cmd:'awk ',file:' -f '}
 keyst = 10B
 printcmd = 'lpr '
 dirlis = {dirlis,cmd:'ls ',sec:'-la ',col:'-1 ',output:'',nohead:'',since:''}
 find = {find,cmd:'find ',ctime:' -ctime ',prune:' -prune ',type:' -type ',  $
         name:' -name ',print:' -print '}
 security = {security,set:'chmod go-rwx ',unset:'chmod a+r,u+w,ug+x,go-w,o-x '}
 quota = 'NA'
;
;  disksep and versep are the ones that are used in DECOMPOSE
;
 syntax = {syntax,redirect:' > ',addon:' >> ',pipe:' | ',background:' & ',   $
           execute:'',allver:'',disksep:':',startpath:'/',midpath:'/',   $
           endpath:'/',versep:';',extsep:'.',listsep:' ',pdelim:' ',   $
           pathlist:':'}
 font = '-ADOBE-Courier-Medium-R-Normal-*-140-75-75-M-90-ISO8859-1'
 syslogin = ''
 trmldev = '/dev/tty'
 nulldev = '/dev/null'
 tmploc = '/tmp'
;
 case (strlowcase(!version.os)) of
    'ultrix':  tapectrl = {tapectrl,device:'/dev/',mt:'mt -f ',back:' bsf ',  $
                           forward:' fsf ',bandf:'',fspace:' fsr ',   $
                           bspace:' bsr ',status:' status',eom:' eom ',   $
                           rewind:' rewind ',weof:' weof 1'}
    'sunos':  tapectrl = {tapectrl,device:'/dev/',mt:'mt -f ',back:' bsf ',  $
                          forward:' fsf ',bandf:' nbsf ',fspace:' fsr ',   $
                          bspace:' bsr ',status:' status',eom:' eom ',   $
                          rewind:' rewind ',weof:' weof 1'}
    'osf':  begin  ; WARNING - these have not been tested
               quota = 'quota '
               tapectrl = {tapectrl,device:'/dev/',mt:'mt -f ',back:' bsf ',  $
                           forward:' fsf ',bandf:'',fspace:' fsr ',   $
                           bspace:' bsr ',status:' status',eom:' soed ',   $
                           rewind:' rewind ',weof:' weof 1'}
            end  ; strlowcase(!version.os) eq 'osf'
;
;  treat all others like SunOS - at least until otherwise determined
;
    else:  tapectrl = {tapectrl,device:'/dev/',mt:'mt -f ',back:' bsf ',  $
                       forward:' fsf ',bandf:' nbsf ',fspace:' fsr ',   $
                       bspace:' bsr ',status:' status',eom:' eom ',  $
                       rewind:' rewind ',weof:' weof 1'}
 endcase  ; strlowcase(!version.os)
;
 tapecopy = {tapecopy,cmd:'dd ',infile:' if=',outfile:' of=',ibs:' ibs=',   $
             obs:' obs=',num:' files=',qic:'512 ',nonqic:'2048 ',   $
             conv:' conv=',sync:'sync '}
;
 return
 end  ; uixplat
