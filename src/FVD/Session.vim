let SessionLoad = 1
if &cp | set nocp | endif
let s:cpo_save=&cpo
set cpo&vim
inoremap <C-Space> 
imap <Nul> <C-Space>
inoremap <expr> <Up> pumvisible() ? "\" : "\<Up>"
inoremap <expr> <S-Tab> pumvisible() ? "\" : "\<S-Tab>"
inoremap <expr> <Down> pumvisible() ? "\" : "\<Down>"
map  h
map <NL> j
map  k
map  l
nnoremap <silent>  :CtrlP
nnoremap   za
nnoremap ,d :YcmShowDetailedDiagnostic
map ,sf z=
map ,sp [s
map ,sn ]s
map ,ss :setlocal spell!
nnoremap ,  :nohlsearch
nnoremap ,f :YcmCompleter FixIt
nnoremap ,g :YcmCompleter GoToImprecise
nnoremap ,s :mksession!
nnoremap ,a :Ag
map ,sa zg
nnoremap ,sv :source $MYVIMRC
nnoremap ,ez :vsp ~/.zshrc
nnoremap ,ev :vsp $MYVIMRC
nnoremap ,u :GundoToggle
nnoremap ,p "0p
map ,cd :cd %:p:h:pwd
map ,te :tabedit =expand("%:p:h")/
nnoremap ,h :split 
nnoremap ,v :vsp 
map ,tk :tabprev
map ,tj :tabnext
map ,ty :tabclose
map ,to :tabonly
map ,tt :tabnew
map ,ba :1,1000 bd!
map ,bd :Bclose
vnoremap // y/" 
vnoremap <silent> 1 :call VisualSelection('f')
vnoremap <silent> 2 :call VisualSelection('b')
nnoremap b ^
nnoremap e $
vmap gx <Plug>NetrwBrowseXVis
nmap gx <Plug>NetrwBrowseX
nnoremap gv `[v`]
nnoremap j gj
nnoremap k gk
vnoremap <silent> <Plug>NetrwBrowseXVis :call netrw#BrowseXVis()
nnoremap <silent> <Plug>NetrwBrowseX :call netrw#BrowseX(expand((exists("g:netrw_gx")? g:netrw_gx : '<cfile>')),netrw#CheckIfRemote())
nmap <F3> :NERDTreeToggle
inoremap <expr> 	 pumvisible() ? "\" : "\	"
vmap ë :m'<-2`>my`<mzgv`yo`z
vmap ê :m'>+`<my`>mzgv`yo`z
nmap ë mz:m-2`z
nmap ê mz:m+`z
inoremap jk 
cabbr q q! 
let &cpo=s:cpo_save
unlet s:cpo_save
set autoread
set backspace=indent,eol,start
set backup
set backupdir=~/.vim-tmp,~/.tmp,~/tmp,/var/tmp,/tmp
set backupskip=/tmp/*,/private/tmp/*
set completefunc=youcompleteme#Complete
set completeopt=preview,menuone
set cpoptions=aAceFsB
set directory=~/.vim-tmp,~/.tmp,~/tmp,/var/tmp,/tmp
set expandtab
set fileencodings=ucs-bom,utf-8,default,latin1
set foldlevelstart=10
set helplang=en
set history=700
set hlsearch
set incsearch
set laststatus=2
set lazyredraw
set mouse=a
set omnifunc=youcompleteme#OmniComplete
set runtimepath=~/.vim,~/.vim/bundle/Vundle.vim,~/.vim/bundle/badwolf,~/.vim/bundle/vim-airline,~/.vim/bundle/vim-airline-themes,~/.vim/bundle/YouCompleteMe,~/.vim/bundle/ctrlp.vim,~/.vim/bundle/nerdtree,~/.vim/bundle/vimtex,/usr/share/vim/vimfiles,/usr/share/vim/vim80/,/usr/share/vim/vimfiles/after,~/.vim/after,~/.vim/bundle/Vundle.vim,~/.vim/bundle/Vundle.vim/after,~/.vim/bundle/badwolf/after,~/.vim/bundle/vim-airline/after,~/.vim/bundle/vim-airline-themes/after,~/.vim/bundle/YouCompleteMe/after,~/.vim/bundle/ctrlp.vim/after,~/.vim/bundle/nerdtree/after,~/.vim/bundle/vimtex/after
set shortmess=filnxtToOc
set showcmd
set showmatch
set showtabline=2
set smartindent
set smarttab
set softtabstop=4
set splitbelow
set splitright
set switchbuf=useopen,usetab,newtab
set tabline=%!airline#extensions#tabline#get()
set tabstop=4
set updatetime=2000
set viminfo=%,'100,<50,s10,h
set wildignore=*.pyc
set wildmenu
let s:so_save = &so | let s:siso_save = &siso | set so=0 siso=0
let v:this_session=expand("<sfile>:p")
silent only
cd ~/pyfy/src/FVD
if expand('%') == '' && !&modified && line('$') <= 1 && getline(1) == ''
  let s:wipebuf = bufnr('%')
endif
set shortmess=aoO
badd +0 Discretize.py
argglobal
silent! argdel *
argadd Discretize.py
edit Discretize.py
set splitbelow splitright
wincmd t
set winheight=1 winwidth=1
argglobal
noremap <buffer>  :s,^\(\s*\)#\s\@!,\1,e:nohlszvj
noremap <buffer> c :s,^\(\s*\)[^# \t]\@=,\1#,e:nohlszvj
setlocal keymap=
setlocal noarabic
setlocal autoindent
setlocal backupcopy=
setlocal balloonexpr=
setlocal nobinary
setlocal nobreakindent
setlocal breakindentopt=
setlocal bufhidden=
setlocal buflisted
setlocal buftype=
setlocal nocindent
setlocal cinkeys=0{,0},0),:,!^F,o,O,e
setlocal cinoptions=
setlocal cinwords=if,else,while,do,for,switch
setlocal colorcolumn=
setlocal comments=b:#,fb:-
setlocal commentstring=#\ %s
setlocal complete=.,w,b,u,t,i
setlocal concealcursor=
setlocal conceallevel=0
setlocal completefunc=youcompleteme#Complete
setlocal nocopyindent
setlocal cryptmethod=
setlocal nocursorbind
setlocal nocursorcolumn
setlocal nocursorline
setlocal define=
setlocal dictionary=
setlocal nodiff
setlocal equalprg=
setlocal errorformat=
setlocal expandtab
if &filetype != 'python'
setlocal filetype=python
endif
setlocal fixendofline
setlocal foldcolumn=0
setlocal foldenable
setlocal foldexpr=0
setlocal foldignore=#
setlocal foldlevel=10
setlocal foldmarker={{{,}}}
set foldmethod=indent
setlocal foldmethod=indent
setlocal foldminlines=1
set foldnestmax=10
setlocal foldnestmax=10
setlocal foldtext=foldtext()
setlocal formatexpr=
setlocal formatoptions=tcq
setlocal formatlistpat=^\\s*\\d\\+[\\]:.)}\\t\ ]\\s*
setlocal grepprg=
setlocal iminsert=2
setlocal imsearch=2
setlocal include=^\\s*\\(from\\|import\\)
setlocal includeexpr=substitute(v:fname,'\\.','/','g')
setlocal indentexpr=GetPythonIndent(v:lnum)
setlocal indentkeys=0{,0},:,!^F,o,O,e,<:>,=elif,=except
setlocal noinfercase
setlocal iskeyword=@,48-57,_,192-255
setlocal keywordprg=pydoc
setlocal nolinebreak
setlocal nolisp
setlocal lispwords=
setlocal nolist
setlocal makeprg=
setlocal matchpairs=(:),{:},[:]
setlocal modeline
setlocal modifiable
setlocal nrformats=bin,octal,hex
set number
setlocal number
setlocal numberwidth=4
setlocal omnifunc=youcompleteme#OmniComplete
setlocal path=
setlocal nopreserveindent
setlocal nopreviewwindow
setlocal quoteescape=\\
setlocal noreadonly
setlocal norelativenumber
setlocal norightleft
setlocal rightleftcmd=search
setlocal noscrollbind
setlocal shiftwidth=4
setlocal noshortname
setlocal signcolumn=auto
setlocal smartindent
setlocal softtabstop=4
setlocal nospell
setlocal spellcapcheck=[.?!]\\_[\\])'\"\	\ ]\\+
setlocal spellfile=
setlocal spelllang=en
setlocal statusline=%!airline#statusline(1)
setlocal suffixesadd=.py
setlocal swapfile
setlocal synmaxcol=3000
if &syntax != 'python'
setlocal syntax=python
endif
setlocal tabstop=8
setlocal tagcase=
setlocal tags=
setlocal textwidth=0
setlocal thesaurus=
setlocal noundofile
setlocal undolevels=-123456
setlocal nowinfixheight
setlocal nowinfixwidth
setlocal wrap
setlocal wrapmargin=0
let s:l = 46 - ((44 * winheight(0) + 25) / 50)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
46
normal! 0
tabnext 1
if exists('s:wipebuf')
  silent exe 'bwipe ' . s:wipebuf
endif
unlet! s:wipebuf
set winheight=1 winwidth=20 shortmess=filnxtToOc
let s:sx = expand("<sfile>:p:r")."x.vim"
if file_readable(s:sx)
  exe "source " . fnameescape(s:sx)
endif
let &so = s:so_save | let &siso = s:siso_save
doautoall SessionLoadPost
unlet SessionLoad
" vim: set ft=vim :
