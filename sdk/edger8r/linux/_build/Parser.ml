type token =
  | EOF
  | TDot
  | TComma
  | TSemicolon
  | TPtr
  | TEqual
  | TLParen
  | TRParen
  | TLBrace
  | TRBrace
  | TLBrack
  | TRBrack
  | Tpublic
  | Tswitchless
  | Tinclude
  | Tconst
  | Tidentifier of (string)
  | Tnumber of (int)
  | Tstring of (string)
  | Tchar
  | Tshort
  | Tunsigned
  | Tint
  | Tfloat
  | Tdouble
  | Tint8
  | Tint16
  | Tint32
  | Tint64
  | Tuint8
  | Tuint16
  | Tuint32
  | Tuint64
  | Tsizet
  | Twchar
  | Tvoid
  | Tlong
  | Tstruct
  | Tunion
  | Tenum
  | Tenclave
  | Tfrom
  | Timport
  | Ttrusted
  | Tuntrusted
  | Tallow
  | Tpropagate_errno

open Parsing;;
let _ = parse_error;;
# 33 "Parser.mly"
open Util				(* for failwithf *)

(* Here we defined some helper routines to check attributes.
 *
 * An alternative approach is to code these rules in Lexer/Parser but
 * it has several drawbacks:
 *
 * 1. Bad extensibility;
 * 2. It grows the table size and down-graded the parsing time;
 * 3. It makes error reporting rigid this way.
 *)

let get_string_from_attr (v: Ast.attr_value) (err_func: int -> string) =
  match v with
      Ast.AString s -> s
    | Ast.ANumber n -> err_func n

(* Check whether 'size' is specified. *)
let has_size (sattr: Ast.ptr_size) =
  sattr.Ast.ps_size <> None
  
(* Check whether 'count' is specified. *)
let has_count (sattr: Ast.ptr_size) =
  sattr.Ast.ps_count <> None

(* Pointers can have the following attributes:
 *
 * 'size'     - specifies the size of the pointer.
 *              e.g. size = 4, size = val ('val' is a parameter);
 *
 * 'count'    - indicates how many of items is managed by the pointer
 *              e.g. count = 100, count = n ('n' is a parameter);
 *
 * 'string'   - indicate the pointer is managing a C string;
 * 'wstring'  - indicate the pointer is managing a wide char string.
 *
 * 'isptr'    - to specify that the foreign type is a pointer.
 * 'isary'    - to specify that the foreign type is an array.
 * 'readonly' - to specify that the foreign type has a 'const' qualifier.
 *
 * 'user_check' - inhibit Edger8r from generating code to check the pointer.
 *
 * 'in'       - the pointer is used as input
 * 'out'      - the pointer is used as output
 *
 * Note that 'size' can be used together with 'count'.
 * 'string' and 'wstring' indicates 'isptr',
 * and they cannot be used with only an 'out' attribute.
 *)
let get_ptr_attr (attr_list: (string * Ast.attr_value) list) =
  let get_new_dir (cds: string) (cda: Ast.ptr_direction) (old: Ast.ptr_direction) =
    if old = Ast.PtrNoDirection then cda
    else if old = Ast.PtrInOut  then failwithf "duplicated attribute: `%s'" cds
    else if old = cda           then failwithf "duplicated attribute: `%s'" cds
    else Ast.PtrInOut
  in
  (* only one 'size' attribute allowed. *)
  let get_new_size (new_value: Ast.attr_value) (old_ptr_size: Ast.ptr_size) =
    if has_size old_ptr_size then
     failwithf "duplicated attribute: `size'"
    else new_value
  in
  (* only one 'count' attribute allowed. *)
  let get_new_count (new_value: Ast.attr_value) (old_ptr_size: Ast.ptr_size) =
    if has_count old_ptr_size then
      failwithf "duplicated attribute: `count'"
    else new_value
  in
  let update_attr (key: string) (value: Ast.attr_value) (res: Ast.ptr_attr) =
    match key with
        "size"     ->
        { res with Ast.pa_size = { res.Ast.pa_size with Ast.ps_size  = Some(get_new_size value res.Ast.pa_size)}}
      | "count"    ->
        { res with Ast.pa_size = { res.Ast.pa_size with Ast.ps_count = Some(get_new_count value res.Ast.pa_size)}}
      | "sizefunc" ->
        failwithf "The attribute 'sizefunc' is deprecated. Please use 'size' attribute instead."
      | "string"  -> { res with Ast.pa_isstr = true; }
      | "wstring" -> { res with Ast.pa_iswstr = true; }
      | "isptr"   -> { res with Ast.pa_isptr = true }
      | "isary"   -> { res with Ast.pa_isary = true }

      | "readonly" -> { res with Ast.pa_rdonly = true }
      | "user_check" -> { res with Ast.pa_chkptr = false }

      | "in"  ->
        let newdir = get_new_dir "in"  Ast.PtrIn  res.Ast.pa_direction
        in { res with Ast.pa_direction = newdir }
      | "out" ->
        let newdir = get_new_dir "out" Ast.PtrOut res.Ast.pa_direction
        in { res with Ast.pa_direction = newdir }
      | _ -> failwithf "unknown attribute: %s" key
  in
  let rec do_get_ptr_attr alist res_attr =
    match alist with
        [] -> res_attr
      | (k,v) :: xs -> do_get_ptr_attr xs (update_attr k v res_attr)
  in
  let has_str_attr (pattr: Ast.ptr_attr) =
    if pattr.Ast.pa_isstr && pattr.Ast.pa_iswstr
    then failwith "`string' and `wstring' are mutual exclusive"
    else (pattr.Ast.pa_isstr || pattr.Ast.pa_iswstr)
  in
  let check_invalid_ptr_size (pattr: Ast.ptr_attr) =
    let ps = pattr.Ast.pa_size in
      if ps <> Ast.empty_ptr_size && has_str_attr pattr
      then failwith "size attributes are mutual exclusive with (w)string attribute"
      else
        if (ps <> Ast.empty_ptr_size || has_str_attr pattr) &&
          pattr.Ast.pa_direction = Ast.PtrNoDirection
        then failwith "size/string attributes must be used with pointer direction"
        else pattr
  in
  let check_ptr_dir (pattr: Ast.ptr_attr) =
    if pattr.Ast.pa_direction <> Ast.PtrNoDirection && pattr.Ast.pa_chkptr = false
    then failwith "pointer direction and `user_check' are mutual exclusive"
    else
      if pattr.Ast.pa_direction = Ast.PtrNoDirection && pattr.Ast.pa_chkptr
      then failwith "pointer/array should have direction attribute or `user_check'"
      else
        if pattr.Ast.pa_direction = Ast.PtrOut && has_str_attr pattr
        then failwith "string/wstring should be used with an `in' attribute"
        else pattr
  in
  let check_invalid_ary_attr (pattr: Ast.ptr_attr) =
    if pattr.Ast.pa_size <> Ast.empty_ptr_size
    then failwith "Pointer size attributes cannot be used with foreign array"
    else
      if not pattr.Ast.pa_isptr
      then
        (* 'pa_chkptr' is default to true unless user specifies 'user_check' *)
        if pattr.Ast.pa_chkptr && pattr.Ast.pa_direction = Ast.PtrNoDirection
        then failwith "array must have direction attribute or `user_check'"
        else pattr
      else
        if has_str_attr pattr
        then failwith "`isary' cannot be used with `string/wstring' together"
        else failwith "`isary' cannot be used with `isptr' together"
  in
  let pattr = do_get_ptr_attr attr_list { Ast.pa_direction = Ast.PtrNoDirection;
                                          Ast.pa_size = Ast.empty_ptr_size;
                                          Ast.pa_isptr = false;
                                          Ast.pa_isary = false;
                                          Ast.pa_isstr = false;
                                          Ast.pa_iswstr = false;
                                          Ast.pa_rdonly = false;
                                          Ast.pa_chkptr = true;
                                        }
  in
    if pattr.Ast.pa_isary
    then check_invalid_ary_attr pattr
    else check_invalid_ptr_size pattr |> check_ptr_dir

(* Untrusted functions can have these attributes:
 *
 * a. 3 mutual exclusive calling convention specifier:
 *     'stdcall', 'fastcall', 'cdecl'.
 *
 * b. 'dllimport' - to import a public symbol.
 *)
let get_func_attr (attr_list: (string * Ast.attr_value) list) =
  let get_new_callconv (key: string) (cur: Ast.call_conv) (old: Ast.call_conv) =
    if old <> Ast.CC_NONE then
      failwithf "unexpected `%s',  conflict with `%s'." key (Ast.get_call_conv_str old)
    else cur
  in
  let update_attr (key: string) (value: Ast.attr_value) (res: Ast.func_attr) =
    match key with
    | "stdcall"  ->
      let callconv = get_new_callconv key Ast.CC_STDCALL res.Ast.fa_convention
      in { res with Ast.fa_convention = callconv}
    | "fastcall" ->
      let callconv = get_new_callconv key Ast.CC_FASTCALL res.Ast.fa_convention
      in { res with Ast.fa_convention = callconv}
    | "cdecl"    ->
      let callconv = get_new_callconv key Ast.CC_CDECL res.Ast.fa_convention
      in { res with Ast.fa_convention = callconv}
    | "dllimport" ->
      if res.Ast.fa_dllimport then failwith "duplicated attribute: `dllimport'"
      else { res with Ast.fa_dllimport = true }
    | _ -> failwithf "invalid function attribute: %s" key
  in
  let rec do_get_func_attr alist res_attr =
    match alist with
      [] -> res_attr
    | (k,v) :: xs -> do_get_func_attr xs (update_attr k v res_attr)
  in do_get_func_attr attr_list { Ast.fa_dllimport = false;
                                  Ast.fa_convention= Ast.CC_NONE;
                                }

(* Some syntax checking against pointer attributes.
 * range: (Lexing.position * Lexing.position)
 *)
let check_ptr_attr (fd: Ast.func_decl) range =
  let fname = fd.Ast.fname in
  let check_const (pattr: Ast.ptr_attr) (identifier: string) =
    let raise_err_direction (direction:string) =
      failwithf "`%s': `%s' is readonly - cannot be used with `%s'"
        fname identifier direction
    in
      if pattr.Ast.pa_rdonly
      then
        match pattr.Ast.pa_direction with
            Ast.PtrOut | Ast.PtrInOut -> raise_err_direction "out"
          | _ -> ()
      else ()
  in
  let check_void_ptr_size (pattr: Ast.ptr_attr) (identifier: string) =
    if pattr.Ast.pa_chkptr && (not (has_size pattr.Ast.pa_size))
    then failwithf "`%s': void pointer `%s' - buffer size unknown" fname identifier
    else ()
  in
  let check_string_ptr_size (atype: Ast.atype) (pattr: Ast.ptr_attr) (identifier: string) =
    if (pattr.Ast.pa_isstr)
    then
      match atype with
      Ast.Ptr(Ast.Char(_)) -> ()
      | _ -> failwithf "`%s': invalid 'string' attribute - `%s' is not char pointer." fname identifier
    else
      if (atype <> Ast.Ptr(Ast.WChar) &&  pattr.Ast.pa_iswstr)
      then failwithf "`%s': invalid 'wstring' attribute - `%s' is not wchar_t pointer." fname identifier
      else ()
  in
  let check_array_dims (atype: Ast.atype) (pattr: Ast.ptr_attr) (declr: Ast.declarator) =
    if Ast.is_array declr then
      if has_size pattr.Ast.pa_size then
        failwithf "`%s': invalid 'size' attribute - `%s' is explicitly declared array." fname declr.Ast.identifier
      else if has_count pattr.Ast.pa_size then
        failwithf "`%s': invalid 'count' attribute - `%s' is explicitly declared array." fname declr.Ast.identifier
      else if pattr.Ast.pa_isary then
        failwithf "`%s': invalid 'isary' attribute - `%s' is explicitly declared array." fname declr.Ast.identifier
    else ()
  in
  let check_pointer_array (atype: Ast.atype) (pattr: Ast.ptr_attr) (declr: Ast.declarator) = 
    let is_ary = (Ast.is_array declr || pattr.Ast.pa_isary) in
    let is_ptr  =
      match atype with
        Ast.Ptr _ -> true
      | _         -> pattr.Ast.pa_isptr
    in
    if is_ary && is_ptr then
      failwithf "`%s': Pointer array not allowed - `%s' is a pointer array." fname declr.Ast.identifier 
    else ()
  in
  let checker (pd: Ast.pdecl) =
    let pt, declr = pd in
    let identifier = declr.Ast.identifier in
      match pt with
          Ast.PTVal _ -> ()
        | Ast.PTPtr(atype, pattr) ->
          if atype = Ast.Ptr(Ast.Void) then (* 'void' pointer, check there is a size or 'user_check' *)
            check_void_ptr_size pattr identifier
          else
            check_pointer_array atype pattr declr;
            check_const pattr identifier;
            check_string_ptr_size atype pattr identifier;
            check_array_dims atype pattr declr
  in
    List.iter checker fd.Ast.plist
# 312 "Parser.ml"
let yytransl_const = [|
    0 (* EOF *);
  257 (* TDot *);
  258 (* TComma *);
  259 (* TSemicolon *);
  260 (* TPtr *);
  261 (* TEqual *);
  262 (* TLParen *);
  263 (* TRParen *);
  264 (* TLBrace *);
  265 (* TRBrace *);
  266 (* TLBrack *);
  267 (* TRBrack *);
  268 (* Tpublic *);
  269 (* Tswitchless *);
  270 (* Tinclude *);
  271 (* Tconst *);
  275 (* Tchar *);
  276 (* Tshort *);
  277 (* Tunsigned *);
  278 (* Tint *);
  279 (* Tfloat *);
  280 (* Tdouble *);
  281 (* Tint8 *);
  282 (* Tint16 *);
  283 (* Tint32 *);
  284 (* Tint64 *);
  285 (* Tuint8 *);
  286 (* Tuint16 *);
  287 (* Tuint32 *);
  288 (* Tuint64 *);
  289 (* Tsizet *);
  290 (* Twchar *);
  291 (* Tvoid *);
  292 (* Tlong *);
  293 (* Tstruct *);
  294 (* Tunion *);
  295 (* Tenum *);
  296 (* Tenclave *);
  297 (* Tfrom *);
  298 (* Timport *);
  299 (* Ttrusted *);
  300 (* Tuntrusted *);
  301 (* Tallow *);
  302 (* Tpropagate_errno *);
    0|]

let yytransl_block = [|
  272 (* Tidentifier *);
  273 (* Tnumber *);
  274 (* Tstring *);
    0|]

let yylhs = "\255\255\
\002\000\002\000\003\000\003\000\004\000\004\000\005\000\005\000\
\006\000\006\000\006\000\006\000\006\000\007\000\007\000\007\000\
\007\000\007\000\007\000\007\000\007\000\007\000\007\000\007\000\
\007\000\007\000\007\000\007\000\007\000\007\000\007\000\007\000\
\007\000\011\000\011\000\012\000\013\000\014\000\014\000\015\000\
\015\000\015\000\016\000\016\000\017\000\017\000\018\000\018\000\
\018\000\018\000\019\000\019\000\020\000\020\000\021\000\021\000\
\021\000\008\000\009\000\010\000\022\000\024\000\025\000\025\000\
\026\000\026\000\027\000\027\000\028\000\028\000\028\000\029\000\
\029\000\029\000\023\000\023\000\030\000\031\000\031\000\032\000\
\033\000\033\000\034\000\035\000\035\000\036\000\036\000\037\000\
\037\000\038\000\038\000\041\000\041\000\042\000\042\000\039\000\
\039\000\040\000\040\000\043\000\043\000\045\000\045\000\045\000\
\046\000\046\000\047\000\048\000\048\000\049\000\049\000\050\000\
\050\000\050\000\044\000\051\000\051\000\051\000\052\000\052\000\
\052\000\052\000\052\000\053\000\001\000\000\000"

let yylen = "\002\000\
\001\000\002\000\001\000\001\000\002\000\003\000\000\000\001\000\
\002\000\003\000\002\000\001\000\001\000\001\000\001\000\001\000\
\001\000\002\000\001\000\001\000\001\000\001\000\001\000\001\000\
\001\000\001\000\001\000\001\000\001\000\001\000\001\000\001\000\
\001\000\001\000\002\000\002\000\003\000\001\000\002\000\001\000\
\001\000\002\000\001\000\002\000\001\000\002\000\002\000\001\000\
\004\000\003\000\002\000\003\000\001\000\003\000\003\000\003\000\
\001\000\002\000\002\000\002\000\004\000\004\000\004\000\004\000\
\000\000\001\000\001\000\003\000\001\000\003\000\003\000\001\000\
\001\000\001\000\002\000\003\000\002\000\001\000\003\000\001\000\
\004\000\004\000\002\000\001\000\002\000\005\000\005\000\001\000\
\002\000\001\000\002\000\000\000\001\000\000\000\001\000\000\000\
\005\000\000\000\003\000\003\000\004\000\002\000\003\000\003\000\
\001\000\003\000\002\000\000\000\001\000\000\000\001\000\000\000\
\002\000\002\000\004\000\000\000\003\000\004\000\000\000\002\000\
\003\000\003\000\002\000\004\000\003\000\002\000"

let yydefred = "\000\000\
\000\000\000\000\000\000\126\000\000\000\119\000\000\000\000\000\
\125\000\124\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\072\000\073\000\074\000\000\000\
\000\000\120\000\123\000\083\000\058\000\059\000\000\000\060\000\
\080\000\000\000\000\000\000\000\000\000\000\000\000\000\122\000\
\121\000\000\000\000\000\000\000\067\000\000\000\084\000\000\000\
\000\000\000\000\000\000\000\000\000\000\033\000\001\000\003\000\
\000\000\016\000\017\000\019\000\020\000\021\000\022\000\023\000\
\024\000\025\000\026\000\027\000\028\000\029\000\000\000\000\000\
\014\000\000\000\012\000\000\000\015\000\000\000\030\000\031\000\
\032\000\000\000\000\000\000\000\000\000\000\000\000\000\063\000\
\000\000\082\000\078\000\000\000\085\000\000\000\000\000\093\000\
\000\000\000\000\000\000\000\000\111\000\000\000\000\000\002\000\
\000\000\008\000\000\000\018\000\005\000\009\000\034\000\000\000\
\000\000\077\000\061\000\000\000\075\000\062\000\064\000\070\000\
\071\000\068\000\000\000\086\000\000\000\000\000\087\000\051\000\
\000\000\000\000\053\000\099\000\000\000\006\000\010\000\035\000\
\000\000\000\000\038\000\000\000\046\000\076\000\079\000\000\000\
\000\000\095\000\000\000\000\000\000\000\052\000\000\000\000\000\
\036\000\000\000\000\000\000\000\039\000\000\000\100\000\000\000\
\097\000\055\000\056\000\054\000\000\000\000\000\000\000\115\000\
\037\000\102\000\000\000\000\000\048\000\000\000\000\000\000\000\
\105\000\101\000\117\000\000\000\109\000\114\000\113\000\000\000\
\103\000\107\000\000\000\047\000\000\000\104\000\118\000\000\000\
\000\000\106\000\000\000"

let yydgoto = "\002\000\
\004\000\073\000\074\000\075\000\076\000\077\000\078\000\079\000\
\080\000\081\000\112\000\138\000\139\000\140\000\141\000\082\000\
\114\000\174\000\101\000\130\000\131\000\021\000\083\000\022\000\
\023\000\043\000\044\000\045\000\024\000\084\000\092\000\034\000\
\025\000\047\000\048\000\027\000\049\000\052\000\050\000\053\000\
\097\000\147\000\126\000\102\000\159\000\176\000\177\000\182\000\
\103\000\168\000\152\000\008\000\005\000"

let yysindex = "\022\000\
\001\255\000\000\062\255\000\000\080\255\000\000\073\000\250\254\
\000\000\000\000\066\255\131\255\164\255\063\255\069\255\173\255\
\196\255\197\255\198\255\199\255\000\000\000\000\000\000\205\255\
\206\255\000\000\000\000\000\000\000\000\000\000\195\255\000\000\
\000\000\170\255\221\255\221\255\142\000\142\000\195\255\000\000\
\000\000\231\255\228\255\236\255\000\000\010\255\000\000\221\255\
\230\255\253\255\221\255\233\255\002\000\000\000\000\000\000\000\
\020\255\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\249\254\224\255\
\000\000\000\000\000\000\245\255\000\000\239\255\000\000\000\000\
\000\000\254\255\032\000\010\000\063\000\006\000\158\255\000\000\
\195\255\000\000\000\000\014\000\000\000\253\255\015\000\000\000\
\142\000\002\000\016\000\100\255\000\000\039\000\142\000\000\000\
\007\000\000\000\023\000\000\000\000\000\000\000\000\000\040\000\
\036\000\000\000\000\000\044\000\000\000\000\000\000\000\000\000\
\000\000\000\000\033\000\000\000\065\255\037\000\000\000\000\000\
\049\000\055\255\000\000\000\000\028\000\000\000\000\000\000\000\
\071\255\064\000\000\000\064\000\000\000\000\000\000\000\069\000\
\060\000\000\000\074\000\160\255\062\000\000\000\075\000\253\254\
\000\000\094\000\068\000\064\000\000\000\070\255\000\000\069\000\
\000\000\000\000\000\000\000\000\051\255\034\000\037\000\000\000\
\000\000\000\000\142\000\097\000\000\000\254\255\118\000\143\255\
\000\000\000\000\000\000\171\255\000\000\000\000\000\000\239\255\
\000\000\000\000\142\000\000\000\093\000\000\000\000\000\040\000\
\239\255\000\000\040\000"

let yyrindex = "\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\098\000\000\000\
\000\000\000\000\101\255\132\255\084\000\084\000\098\000\000\000\
\000\000\067\255\000\000\101\000\000\000\000\000\000\000\101\255\
\000\000\163\255\132\255\000\000\194\255\000\000\000\000\000\000\
\003\255\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\005\255\000\000\
\000\000\008\255\000\000\000\000\000\000\102\255\000\000\000\000\
\000\000\000\000\084\000\000\000\084\000\000\000\000\000\000\000\
\000\000\000\000\000\000\108\000\000\000\225\255\000\000\000\000\
\084\000\001\000\000\000\000\000\000\000\000\000\084\000\000\000\
\005\255\000\000\012\255\000\000\000\000\000\000\000\000\133\255\
\112\255\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\112\000\000\000\000\000\
\061\255\000\000\000\000\000\000\254\254\000\000\000\000\000\000\
\000\000\044\255\000\000\046\255\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\132\000\
\000\000\000\000\000\000\052\255\000\000\084\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\133\000\112\000\000\000\
\000\000\000\000\084\000\048\255\000\000\000\000\084\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\084\000\000\000\084\000\000\000\000\000\124\000\
\000\000\000\000\143\000"

let yygindex = "\000\000\
\000\000\000\000\103\001\000\000\107\001\000\000\089\255\174\001\
\175\001\176\001\137\255\000\000\150\255\047\001\061\001\159\255\
\013\001\000\000\103\255\000\000\039\001\000\000\151\001\000\000\
\000\000\152\001\000\000\101\001\000\000\061\000\027\001\000\000\
\000\000\250\255\157\001\000\000\000\000\000\000\146\001\144\001\
\000\000\029\001\094\001\000\000\038\001\000\000\010\001\000\000\
\000\000\000\000\000\000\000\000\000\000"

let yytablesize = 455
let yytable = "\125\000\
\116\000\026\000\010\000\184\000\175\000\125\000\007\000\011\000\
\004\000\166\000\116\000\013\000\007\000\090\000\004\000\011\000\
\108\000\013\000\007\000\193\000\004\000\011\000\001\000\013\000\
\007\000\091\000\004\000\011\000\109\000\008\000\012\000\013\000\
\014\000\157\000\015\000\175\000\016\000\017\000\104\000\056\000\
\003\000\093\000\167\000\116\000\093\000\041\000\041\000\040\000\
\040\000\157\000\041\000\029\000\040\000\042\000\042\000\105\000\
\149\000\179\000\042\000\041\000\173\000\040\000\057\000\029\000\
\192\000\150\000\091\000\042\000\069\000\006\000\031\000\057\000\
\009\000\195\000\137\000\069\000\170\000\188\000\032\000\100\000\
\144\000\153\000\007\000\028\000\171\000\054\000\033\000\154\000\
\055\000\056\000\057\000\173\000\058\000\059\000\060\000\061\000\
\062\000\063\000\064\000\065\000\066\000\067\000\068\000\069\000\
\172\000\071\000\012\000\013\000\072\000\096\000\128\000\043\000\
\096\000\045\000\045\000\129\000\096\000\043\000\045\000\096\000\
\096\000\096\000\096\000\096\000\096\000\096\000\096\000\096\000\
\096\000\096\000\096\000\096\000\096\000\096\000\096\000\096\000\
\096\000\096\000\096\000\096\000\098\000\098\000\044\000\116\000\
\189\000\116\000\029\000\098\000\044\000\190\000\098\000\098\000\
\098\000\098\000\098\000\098\000\098\000\098\000\098\000\098\000\
\098\000\098\000\098\000\098\000\098\000\098\000\098\000\098\000\
\098\000\098\000\098\000\088\000\123\000\120\000\121\000\162\000\
\163\000\191\000\092\000\030\000\035\000\092\000\092\000\092\000\
\092\000\092\000\092\000\092\000\092\000\092\000\092\000\092\000\
\092\000\092\000\092\000\092\000\092\000\092\000\092\000\092\000\
\092\000\092\000\090\000\036\000\037\000\038\000\039\000\040\000\
\041\000\110\000\042\000\046\000\110\000\110\000\110\000\110\000\
\110\000\110\000\110\000\110\000\110\000\110\000\110\000\110\000\
\110\000\110\000\110\000\110\000\110\000\110\000\110\000\110\000\
\110\000\089\000\011\000\087\000\088\000\089\000\095\000\032\000\
\092\000\099\000\111\000\092\000\092\000\092\000\092\000\092\000\
\092\000\092\000\092\000\092\000\092\000\092\000\092\000\092\000\
\092\000\092\000\092\000\092\000\092\000\092\000\092\000\092\000\
\096\000\091\000\110\000\100\000\117\000\113\000\119\000\123\000\
\110\000\124\000\127\000\110\000\110\000\110\000\110\000\110\000\
\110\000\110\000\110\000\110\000\110\000\110\000\110\000\110\000\
\110\000\110\000\110\000\110\000\110\000\110\000\110\000\110\000\
\115\000\132\000\134\000\136\000\135\000\137\000\142\000\054\000\
\143\000\146\000\055\000\056\000\057\000\148\000\058\000\059\000\
\060\000\061\000\062\000\063\000\064\000\065\000\066\000\067\000\
\068\000\069\000\070\000\071\000\012\000\013\000\072\000\118\000\
\151\000\155\000\158\000\160\000\161\000\129\000\054\000\181\000\
\165\000\055\000\056\000\057\000\154\000\058\000\059\000\060\000\
\061\000\062\000\063\000\064\000\065\000\066\000\067\000\068\000\
\069\000\070\000\071\000\012\000\013\000\072\000\100\000\185\000\
\169\000\007\000\065\000\171\000\054\000\066\000\081\000\055\000\
\056\000\057\000\094\000\058\000\059\000\060\000\061\000\062\000\
\063\000\064\000\065\000\066\000\067\000\068\000\069\000\070\000\
\071\000\012\000\013\000\072\000\187\000\054\000\112\000\108\000\
\055\000\056\000\057\000\050\000\058\000\059\000\060\000\061\000\
\062\000\063\000\064\000\065\000\066\000\067\000\068\000\069\000\
\070\000\071\000\012\000\013\000\072\000\054\000\049\000\106\000\
\055\000\056\000\057\000\107\000\058\000\059\000\060\000\061\000\
\062\000\063\000\064\000\065\000\066\000\067\000\068\000\069\000\
\070\000\071\000\012\000\013\000\072\000\018\000\019\000\020\000\
\156\000\145\000\186\000\164\000\085\000\122\000\086\000\180\000\
\051\000\094\000\098\000\183\000\133\000\178\000\194\000"

let yycheck = "\097\000\
\003\001\008\000\009\001\171\000\158\000\103\000\004\001\014\001\
\004\001\013\001\013\001\004\001\010\001\004\001\010\001\004\001\
\024\001\010\001\016\001\187\000\016\001\010\001\001\000\016\001\
\022\001\016\001\022\001\016\001\036\001\022\001\037\001\038\001\
\039\001\140\000\041\001\189\000\043\001\044\001\019\001\020\001\
\040\001\048\000\046\001\046\001\051\000\002\001\003\001\002\001\
\003\001\156\000\007\001\004\001\007\001\002\001\003\001\036\001\
\002\001\007\001\007\001\016\001\158\000\016\001\002\001\016\001\
\184\000\011\001\016\001\016\001\002\001\008\001\008\001\011\001\
\000\000\193\000\010\001\009\001\007\001\175\000\016\001\010\001\
\016\001\011\001\003\001\018\001\015\001\016\001\018\001\017\001\
\019\001\020\001\021\001\189\000\023\001\024\001\025\001\026\001\
\027\001\028\001\029\001\030\001\031\001\032\001\033\001\034\001\
\035\001\036\001\037\001\038\001\039\001\009\001\011\001\010\001\
\012\001\002\001\003\001\016\001\016\001\016\001\007\001\019\001\
\020\001\021\001\022\001\023\001\024\001\025\001\026\001\027\001\
\028\001\029\001\030\001\031\001\032\001\033\001\034\001\035\001\
\036\001\037\001\038\001\039\001\009\001\010\001\010\001\083\000\
\002\001\085\000\016\001\016\001\016\001\007\001\019\001\020\001\
\021\001\022\001\023\001\024\001\025\001\026\001\027\001\028\001\
\029\001\030\001\031\001\032\001\033\001\034\001\035\001\036\001\
\037\001\038\001\039\001\009\001\002\001\016\001\017\001\016\001\
\017\001\007\001\016\001\016\001\008\001\019\001\020\001\021\001\
\022\001\023\001\024\001\025\001\026\001\027\001\028\001\029\001\
\030\001\031\001\032\001\033\001\034\001\035\001\036\001\037\001\
\038\001\039\001\009\001\008\001\008\001\008\001\008\001\003\001\
\003\001\016\001\016\001\042\001\019\001\020\001\021\001\022\001\
\023\001\024\001\025\001\026\001\027\001\028\001\029\001\030\001\
\031\001\032\001\033\001\034\001\035\001\036\001\037\001\038\001\
\039\001\009\001\014\001\005\001\009\001\002\001\009\001\016\001\
\016\001\009\001\004\001\019\001\020\001\021\001\022\001\023\001\
\024\001\025\001\026\001\027\001\028\001\029\001\030\001\031\001\
\032\001\033\001\034\001\035\001\036\001\037\001\038\001\039\001\
\012\001\009\001\022\001\010\001\003\001\016\001\009\001\002\001\
\016\001\003\001\003\001\019\001\020\001\021\001\022\001\023\001\
\024\001\025\001\026\001\027\001\028\001\029\001\030\001\031\001\
\032\001\033\001\034\001\035\001\036\001\037\001\038\001\039\001\
\009\001\003\001\036\001\004\001\022\001\010\001\003\001\016\001\
\016\001\013\001\019\001\020\001\021\001\005\001\023\001\024\001\
\025\001\026\001\027\001\028\001\029\001\030\001\031\001\032\001\
\033\001\034\001\035\001\036\001\037\001\038\001\039\001\009\001\
\045\001\010\001\006\001\016\001\003\001\016\001\016\001\046\001\
\006\001\019\001\020\001\021\001\017\001\023\001\024\001\025\001\
\026\001\027\001\028\001\029\001\030\001\031\001\032\001\033\001\
\034\001\035\001\036\001\037\001\038\001\039\001\010\001\007\001\
\011\001\022\001\009\001\015\001\016\001\009\001\003\001\019\001\
\020\001\021\001\003\001\023\001\024\001\025\001\026\001\027\001\
\028\001\029\001\030\001\031\001\032\001\033\001\034\001\035\001\
\036\001\037\001\038\001\039\001\015\001\016\001\003\001\003\001\
\019\001\020\001\021\001\016\001\023\001\024\001\025\001\026\001\
\027\001\028\001\029\001\030\001\031\001\032\001\033\001\034\001\
\035\001\036\001\037\001\038\001\039\001\016\001\016\001\057\000\
\019\001\020\001\021\001\057\000\023\001\024\001\025\001\026\001\
\027\001\028\001\029\001\030\001\031\001\032\001\033\001\034\001\
\035\001\036\001\037\001\038\001\039\001\008\000\008\000\008\000\
\138\000\125\000\174\000\149\000\038\000\089\000\039\000\165\000\
\036\000\048\000\051\000\167\000\103\000\160\000\189\000"

let yynames_const = "\
  EOF\000\
  TDot\000\
  TComma\000\
  TSemicolon\000\
  TPtr\000\
  TEqual\000\
  TLParen\000\
  TRParen\000\
  TLBrace\000\
  TRBrace\000\
  TLBrack\000\
  TRBrack\000\
  Tpublic\000\
  Tswitchless\000\
  Tinclude\000\
  Tconst\000\
  Tchar\000\
  Tshort\000\
  Tunsigned\000\
  Tint\000\
  Tfloat\000\
  Tdouble\000\
  Tint8\000\
  Tint16\000\
  Tint32\000\
  Tint64\000\
  Tuint8\000\
  Tuint16\000\
  Tuint32\000\
  Tuint64\000\
  Tsizet\000\
  Twchar\000\
  Tvoid\000\
  Tlong\000\
  Tstruct\000\
  Tunion\000\
  Tenum\000\
  Tenclave\000\
  Tfrom\000\
  Timport\000\
  Ttrusted\000\
  Tuntrusted\000\
  Tallow\000\
  Tpropagate_errno\000\
  "

let yynames_block = "\
  Tidentifier\000\
  Tnumber\000\
  Tstring\000\
  "

let yyact = [|
  (fun _ -> failwith "parser")
; (fun __caml_parser_env ->
    Obj.repr(
# 320 "Parser.mly"
                 ( Ast.Char Ast.Signed )
# 679 "Parser.ml"
               : 'char_type))
; (fun __caml_parser_env ->
    Obj.repr(
# 321 "Parser.mly"
                    ( Ast.Char Ast.Unsigned )
# 685 "Parser.ml"
               : 'char_type))
; (fun __caml_parser_env ->
    Obj.repr(
# 325 "Parser.mly"
                     ( Ast.IShort )
# 691 "Parser.ml"
               : 'ex_shortness))
; (fun __caml_parser_env ->
    Obj.repr(
# 326 "Parser.mly"
          ( Ast.ILong )
# 697 "Parser.ml"
               : 'ex_shortness))
; (fun __caml_parser_env ->
    Obj.repr(
# 329 "Parser.mly"
                          ( Ast.LLong Ast.Signed )
# 703 "Parser.ml"
               : 'longlong))
; (fun __caml_parser_env ->
    Obj.repr(
# 330 "Parser.mly"
                          ( Ast.LLong Ast.Unsigned )
# 709 "Parser.ml"
               : 'longlong))
; (fun __caml_parser_env ->
    Obj.repr(
# 332 "Parser.mly"
                       ( Ast.INone )
# 715 "Parser.ml"
               : 'shortness))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'ex_shortness) in
    Obj.repr(
# 333 "Parser.mly"
                 ( _1 )
# 722 "Parser.ml"
               : 'shortness))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'shortness) in
    Obj.repr(
# 336 "Parser.mly"
                         (
      Ast.Int { Ast.ia_signedness = Ast.Signed; Ast.ia_shortness = _1 }
    )
# 731 "Parser.ml"
               : 'int_type))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'shortness) in
    Obj.repr(
# 339 "Parser.mly"
                             (
      Ast.Int { Ast.ia_signedness = Ast.Unsigned; Ast.ia_shortness = _2 }
    )
# 740 "Parser.ml"
               : 'int_type))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'shortness) in
    Obj.repr(
# 342 "Parser.mly"
                        (
      Ast.Int { Ast.ia_signedness = Ast.Unsigned; Ast.ia_shortness = _2 }
    )
# 749 "Parser.ml"
               : 'int_type))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'longlong) in
    Obj.repr(
# 345 "Parser.mly"
             ( _1 )
# 756 "Parser.ml"
               : 'int_type))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'ex_shortness) in
    Obj.repr(
# 346 "Parser.mly"
                 (
      Ast.Int { Ast.ia_signedness = Ast.Signed; Ast.ia_shortness = _1 }
    )
# 765 "Parser.ml"
               : 'int_type))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'char_type) in
    Obj.repr(
# 352 "Parser.mly"
              ( _1 )
# 772 "Parser.ml"
               : 'type_spec))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'int_type) in
    Obj.repr(
# 353 "Parser.mly"
              ( _1 )
# 779 "Parser.ml"
               : 'type_spec))
; (fun __caml_parser_env ->
    Obj.repr(
# 355 "Parser.mly"
             ( Ast.Float )
# 785 "Parser.ml"
               : 'type_spec))
; (fun __caml_parser_env ->
    Obj.repr(
# 356 "Parser.mly"
             ( Ast.Double )
# 791 "Parser.ml"
               : 'type_spec))
; (fun __caml_parser_env ->
    Obj.repr(
# 357 "Parser.mly"
                  ( Ast.LDouble )
# 797 "Parser.ml"
               : 'type_spec))
; (fun __caml_parser_env ->
    Obj.repr(
# 359 "Parser.mly"
             ( Ast.Int8 )
# 803 "Parser.ml"
               : 'type_spec))
; (fun __caml_parser_env ->
    Obj.repr(
# 360 "Parser.mly"
             ( Ast.Int16 )
# 809 "Parser.ml"
               : 'type_spec))
; (fun __caml_parser_env ->
    Obj.repr(
# 361 "Parser.mly"
             ( Ast.Int32 )
# 815 "Parser.ml"
               : 'type_spec))
; (fun __caml_parser_env ->
    Obj.repr(
# 362 "Parser.mly"
             ( Ast.Int64 )
# 821 "Parser.ml"
               : 'type_spec))
; (fun __caml_parser_env ->
    Obj.repr(
# 363 "Parser.mly"
             ( Ast.UInt8 )
# 827 "Parser.ml"
               : 'type_spec))
; (fun __caml_parser_env ->
    Obj.repr(
# 364 "Parser.mly"
             ( Ast.UInt16 )
# 833 "Parser.ml"
               : 'type_spec))
; (fun __caml_parser_env ->
    Obj.repr(
# 365 "Parser.mly"
             ( Ast.UInt32 )
# 839 "Parser.ml"
               : 'type_spec))
; (fun __caml_parser_env ->
    Obj.repr(
# 366 "Parser.mly"
             ( Ast.UInt64 )
# 845 "Parser.ml"
               : 'type_spec))
; (fun __caml_parser_env ->
    Obj.repr(
# 367 "Parser.mly"
             ( Ast.SizeT )
# 851 "Parser.ml"
               : 'type_spec))
; (fun __caml_parser_env ->
    Obj.repr(
# 368 "Parser.mly"
             ( Ast.WChar )
# 857 "Parser.ml"
               : 'type_spec))
; (fun __caml_parser_env ->
    Obj.repr(
# 369 "Parser.mly"
             ( Ast.Void )
# 863 "Parser.ml"
               : 'type_spec))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'struct_specifier) in
    Obj.repr(
# 371 "Parser.mly"
                     ( _1 )
# 870 "Parser.ml"
               : 'type_spec))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'union_specifier) in
    Obj.repr(
# 372 "Parser.mly"
                     ( _1 )
# 877 "Parser.ml"
               : 'type_spec))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'enum_specifier) in
    Obj.repr(
# 373 "Parser.mly"
                     ( _1 )
# 884 "Parser.ml"
               : 'type_spec))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : string) in
    Obj.repr(
# 374 "Parser.mly"
                     ( Ast.Foreign(_1) )
# 891 "Parser.ml"
               : 'type_spec))
; (fun __caml_parser_env ->
    Obj.repr(
# 377 "Parser.mly"
                 ( fun ii -> Ast.Ptr(ii) )
# 897 "Parser.ml"
               : 'pointer))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'pointer) in
    Obj.repr(
# 378 "Parser.mly"
                 ( fun ii -> Ast.Ptr(_1 ii) )
# 904 "Parser.ml"
               : 'pointer))
; (fun __caml_parser_env ->
    Obj.repr(
# 381 "Parser.mly"
                                         ( failwith "Flexible array is not supported." )
# 910 "Parser.ml"
               : 'empty_dimension))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : int) in
    Obj.repr(
# 382 "Parser.mly"
                                         ( if _2 <> 0 then [_2]
                                           else failwith "Zero-length array is not supported." )
# 918 "Parser.ml"
               : 'fixed_dimension))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'fixed_dimension) in
    Obj.repr(
# 385 "Parser.mly"
                                     ( _1 )
# 925 "Parser.ml"
               : 'fixed_size_array))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'fixed_size_array) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'fixed_dimension) in
    Obj.repr(
# 386 "Parser.mly"
                                     ( _1 @ _2 )
# 933 "Parser.ml"
               : 'fixed_size_array))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'fixed_size_array) in
    Obj.repr(
# 389 "Parser.mly"
                                     ( _1 )
# 940 "Parser.ml"
               : 'array_size))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'empty_dimension) in
    Obj.repr(
# 390 "Parser.mly"
                                     ( _1 )
# 947 "Parser.ml"
               : 'array_size))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'empty_dimension) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'fixed_size_array) in
    Obj.repr(
# 391 "Parser.mly"
                                     ( _1 @ _2 )
# 955 "Parser.ml"
               : 'array_size))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'type_spec) in
    Obj.repr(
# 394 "Parser.mly"
                      ( _1 )
# 962 "Parser.ml"
               : 'all_type))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'type_spec) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'pointer) in
    Obj.repr(
# 395 "Parser.mly"
                      ( _2 _1 )
# 970 "Parser.ml"
               : 'all_type))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : string) in
    Obj.repr(
# 398 "Parser.mly"
                           ( { Ast.identifier = _1; Ast.array_dims = []; } )
# 977 "Parser.ml"
               : 'declarator))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : string) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'array_size) in
    Obj.repr(
# 399 "Parser.mly"
                           ( { Ast.identifier = _1; Ast.array_dims = _2; } )
# 985 "Parser.ml"
               : 'declarator))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'attr_block) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'all_type) in
    Obj.repr(
# 408 "Parser.mly"
                                (
    let attr = get_ptr_attr _1 in
    (*check the type is build in type or used defined type.*)
    let rec is_foreign s =
      match s with
        Ast.Ptr(a) -> is_foreign a
      | Ast.Foreign _ -> true
      | _ -> false
    in
    let is_bare_foreign s =
      match s with
      | Ast.Foreign _ -> true
      | _ -> false
    in
    (*'isptr', 'isary', only allowed for bare user defined type.*)
    (*'readonly' only allowed for user defined type.*)
    if attr.Ast.pa_isptr && not (is_bare_foreign _2) then
      failwithf "'isptr', attributes are only for user defined type, not for `%s'." (Ast.get_tystr _2)
    else if attr.Ast.pa_isary && not (is_bare_foreign _2) then
      failwithf "'isary', attributes are only for user defined type, not for `%s'." (Ast.get_tystr _2)
    else if attr.Ast.pa_rdonly && not (is_foreign _2) then
      failwithf "'readonly', attributes are only for user defined type, not for `%s'." (Ast.get_tystr _2)
    else if attr.Ast.pa_rdonly && not (attr.Ast.pa_isptr) then
      failwithf "'readonly' attribute is only used with 'isptr' attribute."    else
    match _2 with
      Ast.Ptr _ -> fun x -> Ast.PTPtr(_2, get_ptr_attr _1)
    | _         ->
      if _1 <> [] then
        let attr = get_ptr_attr _1 in
        match _2 with
          Ast.Foreign s ->
            if attr.Ast.pa_isptr || attr.Ast.pa_isary then fun x -> Ast.PTPtr(_2, attr)
            else
              (* thinking about 'user_defined_type var[4]' *)
              fun is_ary ->
                if is_ary then Ast.PTPtr(_2, attr)
                else failwithf "`%s' is considered plain type but decorated with pointer attributes" s
        | _ ->
          fun is_ary ->
            if is_ary then Ast.PTPtr(_2, attr)
            else failwithf "unexpected pointer attributes for `%s'" (Ast.get_tystr _2)
      else
        fun is_ary ->
          if is_ary then Ast.PTPtr(_2, get_ptr_attr [])
          else  Ast.PTVal _2
    )
# 1038 "Parser.ml"
               : 'param_type))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'all_type) in
    Obj.repr(
# 454 "Parser.mly"
             (
    match _1 with
      Ast.Ptr _ -> fun x -> Ast.PTPtr(_1, get_ptr_attr [])
    | _         ->
      fun is_ary ->
        if is_ary then Ast.PTPtr(_1, get_ptr_attr [])
        else  Ast.PTVal _1
    )
# 1052 "Parser.ml"
               : 'param_type))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 3 : 'attr_block) in
    let _3 = (Parsing.peek_val __caml_parser_env 1 : 'type_spec) in
    let _4 = (Parsing.peek_val __caml_parser_env 0 : 'pointer) in
    Obj.repr(
# 462 "Parser.mly"
                                        (
      let attr = get_ptr_attr _1
      in fun x -> Ast.PTPtr(_4 _3, { attr with Ast.pa_rdonly = true })
    )
# 1064 "Parser.ml"
               : 'param_type))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'type_spec) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'pointer) in
    Obj.repr(
# 466 "Parser.mly"
                             (
      let attr = get_ptr_attr []
      in fun x -> Ast.PTPtr(_3 _2, { attr with Ast.pa_rdonly = true })
    )
# 1075 "Parser.ml"
               : 'param_type))
; (fun __caml_parser_env ->
    Obj.repr(
# 473 "Parser.mly"
                                  ( failwith "no attribute specified." )
# 1081 "Parser.ml"
               : 'attr_block))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'key_val_pairs) in
    Obj.repr(
# 474 "Parser.mly"
                                  ( _2 )
# 1088 "Parser.ml"
               : 'attr_block))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'key_val_pair) in
    Obj.repr(
# 477 "Parser.mly"
                                      ( [_1] )
# 1095 "Parser.ml"
               : 'key_val_pairs))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'key_val_pairs) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'key_val_pair) in
    Obj.repr(
# 478 "Parser.mly"
                                      (  _3 :: _1 )
# 1103 "Parser.ml"
               : 'key_val_pairs))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : string) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : string) in
    Obj.repr(
# 481 "Parser.mly"
                                             ( (_1, Ast.AString(_3)) )
# 1111 "Parser.ml"
               : 'key_val_pair))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : string) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : int) in
    Obj.repr(
# 482 "Parser.mly"
                                             ( (_1, Ast.ANumber(_3)) )
# 1119 "Parser.ml"
               : 'key_val_pair))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : string) in
    Obj.repr(
# 483 "Parser.mly"
                                             ( (_1, Ast.AString("")) )
# 1126 "Parser.ml"
               : 'key_val_pair))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 0 : string) in
    Obj.repr(
# 486 "Parser.mly"
                                      ( Ast.Struct(_2) )
# 1133 "Parser.ml"
               : 'struct_specifier))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 0 : string) in
    Obj.repr(
# 487 "Parser.mly"
                                      ( Ast.Union(_2) )
# 1140 "Parser.ml"
               : 'union_specifier))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 0 : string) in
    Obj.repr(
# 488 "Parser.mly"
                                      ( Ast.Enum(_2) )
# 1147 "Parser.ml"
               : 'enum_specifier))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 3 : 'struct_specifier) in
    let _3 = (Parsing.peek_val __caml_parser_env 1 : 'member_list) in
    Obj.repr(
# 490 "Parser.mly"
                                                                (
    let s = { Ast.sname = (match _1 with Ast.Struct s -> s | _ -> "");
              Ast.mlist = List.rev _3; }
    in Ast.StructDef(s)
  )
# 1159 "Parser.ml"
               : 'struct_definition))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 3 : 'union_specifier) in
    let _3 = (Parsing.peek_val __caml_parser_env 1 : 'member_list) in
    Obj.repr(
# 496 "Parser.mly"
                                                              (
    let s = { Ast.sname = (match _1 with Ast.Union s -> s | _ -> "");
              Ast.mlist = List.rev _3; }
    in Ast.UnionDef(s)
  )
# 1171 "Parser.ml"
               : 'union_definition))
; (fun __caml_parser_env ->
    let _3 = (Parsing.peek_val __caml_parser_env 1 : 'enum_body) in
    Obj.repr(
# 503 "Parser.mly"
                                                 (
      let e = { Ast.enname = ""; Ast.enbody = _3; }
      in Ast.EnumDef(e)
    )
# 1181 "Parser.ml"
               : 'enum_definition))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 3 : 'enum_specifier) in
    let _3 = (Parsing.peek_val __caml_parser_env 1 : 'enum_body) in
    Obj.repr(
# 507 "Parser.mly"
                                             (
      let e = { Ast.enname = (match _1 with Ast.Enum s -> s | _ -> "");
                Ast.enbody = _3; }
      in Ast.EnumDef(e)
    )
# 1193 "Parser.ml"
               : 'enum_definition))
; (fun __caml_parser_env ->
    Obj.repr(
# 514 "Parser.mly"
                       ( [] )
# 1199 "Parser.ml"
               : 'enum_body))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'enum_eles) in
    Obj.repr(
# 515 "Parser.mly"
                       ( List.rev _1 )
# 1206 "Parser.ml"
               : 'enum_body))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'enum_ele) in
    Obj.repr(
# 518 "Parser.mly"
                              ( [_1] )
# 1213 "Parser.ml"
               : 'enum_eles))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'enum_eles) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'enum_ele) in
    Obj.repr(
# 519 "Parser.mly"
                              ( _3 :: _1 )
# 1221 "Parser.ml"
               : 'enum_eles))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : string) in
    Obj.repr(
# 522 "Parser.mly"
                                   ( (_1, Ast.EnumValNone) )
# 1228 "Parser.ml"
               : 'enum_ele))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : string) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : string) in
    Obj.repr(
# 523 "Parser.mly"
                                   ( (_1, Ast.EnumVal (Ast.AString _3)) )
# 1236 "Parser.ml"
               : 'enum_ele))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : string) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : int) in
    Obj.repr(
# 524 "Parser.mly"
                                   ( (_1, Ast.EnumVal (Ast.ANumber _3)) )
# 1244 "Parser.ml"
               : 'enum_ele))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'struct_definition) in
    Obj.repr(
# 527 "Parser.mly"
                                      ( _1 )
# 1251 "Parser.ml"
               : 'composite_defs))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'union_definition) in
    Obj.repr(
# 528 "Parser.mly"
                                      ( _1 )
# 1258 "Parser.ml"
               : 'composite_defs))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'enum_definition) in
    Obj.repr(
# 529 "Parser.mly"
                                      ( _1 )
# 1265 "Parser.ml"
               : 'composite_defs))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'member_def) in
    Obj.repr(
# 532 "Parser.mly"
                                      ( [_1] )
# 1272 "Parser.ml"
               : 'member_list))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'member_list) in
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'member_def) in
    Obj.repr(
# 533 "Parser.mly"
                                      ( _2 :: _1 )
# 1280 "Parser.ml"
               : 'member_list))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'all_type) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'declarator) in
    Obj.repr(
# 536 "Parser.mly"
                                ( (_1, _2) )
# 1288 "Parser.ml"
               : 'member_def))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : string) in
    Obj.repr(
# 541 "Parser.mly"
                                  ( [_1] )
# 1295 "Parser.ml"
               : 'func_list))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'func_list) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : string) in
    Obj.repr(
# 542 "Parser.mly"
                                  ( _3 :: _1 )
# 1303 "Parser.ml"
               : 'func_list))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : string) in
    Obj.repr(
# 545 "Parser.mly"
                                  ( _1 )
# 1310 "Parser.ml"
               : 'module_path))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 2 : 'module_path) in
    let _4 = (Parsing.peek_val __caml_parser_env 0 : 'func_list) in
    Obj.repr(
# 547 "Parser.mly"
                                                         (
      { Ast.mname = _2; Ast.flist = List.rev _4; }
    )
# 1320 "Parser.ml"
               : 'import_declaration))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 2 : 'module_path) in
    Obj.repr(
# 550 "Parser.mly"
                                   (
      { Ast.mname = _2; Ast.flist = ["*"]; }
    )
# 1329 "Parser.ml"
               : 'import_declaration))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 0 : string) in
    Obj.repr(
# 555 "Parser.mly"
                                      ( _2 )
# 1336 "Parser.ml"
               : 'include_declaration))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'include_declaration) in
    Obj.repr(
# 557 "Parser.mly"
                                             ( [_1] )
# 1343 "Parser.ml"
               : 'include_declarations))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'include_declarations) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'include_declaration) in
    Obj.repr(
# 558 "Parser.mly"
                                             ( _2 :: _1 )
# 1351 "Parser.ml"
               : 'include_declarations))
; (fun __caml_parser_env ->
    let _3 = (Parsing.peek_val __caml_parser_env 2 : 'trusted_block) in
    Obj.repr(
# 564 "Parser.mly"
                                                                     (
      List.rev _3
    )
# 1360 "Parser.ml"
               : 'enclave_functions))
; (fun __caml_parser_env ->
    let _3 = (Parsing.peek_val __caml_parser_env 2 : 'untrusted_block) in
    Obj.repr(
# 567 "Parser.mly"
                                                          (
      List.rev _3
    )
# 1369 "Parser.ml"
               : 'enclave_functions))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'trusted_functions) in
    Obj.repr(
# 572 "Parser.mly"
                                             ( _1 )
# 1376 "Parser.ml"
               : 'trusted_block))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'include_declarations) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'trusted_functions) in
    Obj.repr(
# 573 "Parser.mly"
                                             (
      trusted_headers := !trusted_headers @ List.rev _1; _2
    )
# 1386 "Parser.ml"
               : 'trusted_block))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'untrusted_functions) in
    Obj.repr(
# 578 "Parser.mly"
                                             ( _1 )
# 1393 "Parser.ml"
               : 'untrusted_block))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'include_declarations) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'untrusted_functions) in
    Obj.repr(
# 579 "Parser.mly"
                                             (
      untrusted_headers := !untrusted_headers @ List.rev _1; _2
    )
# 1403 "Parser.ml"
               : 'untrusted_block))
; (fun __caml_parser_env ->
    Obj.repr(
# 585 "Parser.mly"
                               ( true )
# 1409 "Parser.ml"
               : 'access_modifier))
; (fun __caml_parser_env ->
    Obj.repr(
# 586 "Parser.mly"
                               ( false  )
# 1415 "Parser.ml"
               : 'access_modifier))
; (fun __caml_parser_env ->
    Obj.repr(
# 590 "Parser.mly"
                                     ( false )
# 1421 "Parser.ml"
               : 'switchless_annotation))
; (fun __caml_parser_env ->
    Obj.repr(
# 591 "Parser.mly"
                                     ( true  )
# 1427 "Parser.ml"
               : 'switchless_annotation))
; (fun __caml_parser_env ->
    Obj.repr(
# 594 "Parser.mly"
                                          ( [] )
# 1433 "Parser.ml"
               : 'trusted_functions))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 4 : 'trusted_functions) in
    let _2 = (Parsing.peek_val __caml_parser_env 3 : 'access_modifier) in
    let _3 = (Parsing.peek_val __caml_parser_env 2 : 'func_def) in
    let _4 = (Parsing.peek_val __caml_parser_env 1 : 'switchless_annotation) in
    Obj.repr(
# 595 "Parser.mly"
                                                                                (
      check_ptr_attr _3 (symbol_start_pos(), symbol_end_pos());
      Ast.Trusted { Ast.tf_fdecl = _3; Ast.tf_is_priv = _2; Ast.tf_is_switchless = _4 } :: _1
    )
# 1446 "Parser.ml"
               : 'trusted_functions))
; (fun __caml_parser_env ->
    Obj.repr(
# 601 "Parser.mly"
                                                      ( [] )
# 1452 "Parser.ml"
               : 'untrusted_functions))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'untrusted_functions) in
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'untrusted_func_def) in
    Obj.repr(
# 602 "Parser.mly"
                                                      ( _2 :: _1 )
# 1460 "Parser.ml"
               : 'untrusted_functions))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'all_type) in
    let _2 = (Parsing.peek_val __caml_parser_env 1 : string) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'parameter_list) in
    Obj.repr(
# 605 "Parser.mly"
                                              (
      { Ast.fname = _2; Ast.rtype = _1; Ast.plist = List.rev _3 ; }
    )
# 1471 "Parser.ml"
               : 'func_def))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 3 : 'all_type) in
    let _2 = (Parsing.peek_val __caml_parser_env 2 : 'array_size) in
    let _3 = (Parsing.peek_val __caml_parser_env 1 : string) in
    let _4 = (Parsing.peek_val __caml_parser_env 0 : 'parameter_list) in
    Obj.repr(
# 608 "Parser.mly"
                                                   (
      failwithf "%s: returning an array is not supported - use pointer instead." _3
    )
# 1483 "Parser.ml"
               : 'func_def))
; (fun __caml_parser_env ->
    Obj.repr(
# 613 "Parser.mly"
                                   ( [] )
# 1489 "Parser.ml"
               : 'parameter_list))
; (fun __caml_parser_env ->
    Obj.repr(
# 614 "Parser.mly"
                                   ( [] )
# 1495 "Parser.ml"
               : 'parameter_list))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'parameter_defs) in
    Obj.repr(
# 615 "Parser.mly"
                                   ( _2 )
# 1502 "Parser.ml"
               : 'parameter_list))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'parameter_def) in
    Obj.repr(
# 618 "Parser.mly"
                                        ( [_1] )
# 1509 "Parser.ml"
               : 'parameter_defs))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'parameter_defs) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'parameter_def) in
    Obj.repr(
# 619 "Parser.mly"
                                        ( _3 :: _1 )
# 1517 "Parser.ml"
               : 'parameter_defs))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'param_type) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'declarator) in
    Obj.repr(
# 622 "Parser.mly"
                                     (
    let pt = _1 (Ast.is_array _2) in
    let is_void =
      match pt with
          Ast.PTVal v -> v = Ast.Void
        | _           -> false
    in
      if is_void then
        failwithf "parameter `%s' has `void' type." _2.Ast.identifier
      else
        (pt, _2)
  )
# 1536 "Parser.ml"
               : 'parameter_def))
; (fun __caml_parser_env ->
    Obj.repr(
# 636 "Parser.mly"
                               ( false )
# 1542 "Parser.ml"
               : 'propagate_errno))
; (fun __caml_parser_env ->
    Obj.repr(
# 637 "Parser.mly"
                               ( true  )
# 1548 "Parser.ml"
               : 'propagate_errno))
; (fun __caml_parser_env ->
    Obj.repr(
# 640 "Parser.mly"
                                  ( [] )
# 1554 "Parser.ml"
               : 'untrusted_prefixes))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'attr_block) in
    Obj.repr(
# 641 "Parser.mly"
                         ( _1  )
# 1561 "Parser.ml"
               : 'untrusted_prefixes))
; (fun __caml_parser_env ->
    Obj.repr(
# 644 "Parser.mly"
                                     (  (false, false) )
# 1567 "Parser.ml"
               : 'untrusted_postfixes))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'switchless_annotation) in
    Obj.repr(
# 645 "Parser.mly"
                                            ( (true, _2) )
# 1574 "Parser.ml"
               : 'untrusted_postfixes))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'propagate_errno) in
    Obj.repr(
# 646 "Parser.mly"
                                 ( (_2, true) )
# 1581 "Parser.ml"
               : 'untrusted_postfixes))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 3 : 'untrusted_prefixes) in
    let _2 = (Parsing.peek_val __caml_parser_env 2 : 'func_def) in
    let _3 = (Parsing.peek_val __caml_parser_env 1 : 'allow_list) in
    let _4 = (Parsing.peek_val __caml_parser_env 0 : 'untrusted_postfixes) in
    Obj.repr(
# 649 "Parser.mly"
                                                                               (
      check_ptr_attr _2 (symbol_start_pos(), symbol_end_pos());
      let fattr = get_func_attr _1 in
      Ast.Untrusted { Ast.uf_fdecl = _2; Ast.uf_fattr = fattr; Ast.uf_allow_list = _3; Ast.uf_propagate_errno = fst _4; Ast.uf_is_switchless = snd _4; }
    )
# 1595 "Parser.ml"
               : 'untrusted_func_def))
; (fun __caml_parser_env ->
    Obj.repr(
# 656 "Parser.mly"
                                     ( [] )
# 1601 "Parser.ml"
               : 'allow_list))
; (fun __caml_parser_env ->
    Obj.repr(
# 657 "Parser.mly"
                                     ( [] )
# 1607 "Parser.ml"
               : 'allow_list))
; (fun __caml_parser_env ->
    let _3 = (Parsing.peek_val __caml_parser_env 1 : 'func_list) in
    Obj.repr(
# 658 "Parser.mly"
                                     ( _3 )
# 1614 "Parser.ml"
               : 'allow_list))
; (fun __caml_parser_env ->
    Obj.repr(
# 664 "Parser.mly"
                           ( [] )
# 1620 "Parser.ml"
               : 'expressions))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'expressions) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'include_declaration) in
    Obj.repr(
# 665 "Parser.mly"
                                              ( Ast.Include(_2)   :: _1 )
# 1628 "Parser.ml"
               : 'expressions))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'expressions) in
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'import_declaration) in
    Obj.repr(
# 666 "Parser.mly"
                                              ( Ast.Importing(_2) :: _1 )
# 1636 "Parser.ml"
               : 'expressions))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'expressions) in
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'composite_defs) in
    Obj.repr(
# 667 "Parser.mly"
                                              ( Ast.Composite(_2) :: _1 )
# 1644 "Parser.ml"
               : 'expressions))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'expressions) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'enclave_functions) in
    Obj.repr(
# 668 "Parser.mly"
                                              ( Ast.Interface(_2) :: _1 )
# 1652 "Parser.ml"
               : 'expressions))
; (fun __caml_parser_env ->
    let _3 = (Parsing.peek_val __caml_parser_env 1 : 'expressions) in
    Obj.repr(
# 671 "Parser.mly"
                                                  (
      { Ast.ename = "";
        Ast.eexpr = List.rev _3 }
    )
# 1662 "Parser.ml"
               : 'enclave_def))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'enclave_def) in
    Obj.repr(
# 680 "Parser.mly"
                                          ( _1 )
# 1669 "Parser.ml"
               : Ast.enclave))
(* Entry start_parsing *)
; (fun __caml_parser_env -> raise (Parsing.YYexit (Parsing.peek_val __caml_parser_env 0)))
|]
let yytables =
  { Parsing.actions=yyact;
    Parsing.transl_const=yytransl_const;
    Parsing.transl_block=yytransl_block;
    Parsing.lhs=yylhs;
    Parsing.len=yylen;
    Parsing.defred=yydefred;
    Parsing.dgoto=yydgoto;
    Parsing.sindex=yysindex;
    Parsing.rindex=yyrindex;
    Parsing.gindex=yygindex;
    Parsing.tablesize=yytablesize;
    Parsing.table=yytable;
    Parsing.check=yycheck;
    Parsing.error_function=parse_error;
    Parsing.names_const=yynames_const;
    Parsing.names_block=yynames_block }
let start_parsing (lexfun : Lexing.lexbuf -> token) (lexbuf : Lexing.lexbuf) =
   (Parsing.yyparse yytables 1 lexfun lexbuf : Ast.enclave)
;;
