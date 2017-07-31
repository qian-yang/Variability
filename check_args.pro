;+
; NAME:
;   check_args
; PURPOSE:
;   check whether a tag name exist in a structure, and whether it is ture or not
;
; INPUT:
;   struct: structure
;   name: tag name
;
; OPTIONAL INPUTS:
;
; OUTPUT:
;   res: 1 if true, 0 if not ture
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   depend on the TAG_EXIST function
;
; REVISION HISTORY:
;   29-Jul-2017  Written by Qian Yang, qianyang.astro@gmail.com
;
; Tudo:
;-

function check_args, struct, name
  res = 0
  if tag_exist(struct, name, index=id) then begin
    if struct.(id)[0] then res = 1
  endif
  return, res
end
