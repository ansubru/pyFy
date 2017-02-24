"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" Maintainer: 
" "       ask
" "       ananda@chalmers.se
" "
" " Version: 
" "       1.0 - 19/12/16 12:00:00
" "
" " Sections:
" "    -> General
" "    -> Colors and Fonts
" "    -> Tabs and spacing
" "    -> UI Config
" "    -> Visual mode related
" "    -> Searching
" "    -> Folding
" "    -> Movement, tabs, windows and buffers
" "    -> Leader shortcuts
" "    -> CtrlP settings
" "    -> Launch config
" "    -> Backups
" "    -> Spell check
" "    -> Autogroup
" "
" " 
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" " => General
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
set nocompatible              " be iMproved, required
filetype off                  " required

"set the runtime path to include Vundle and initialize
set rtp+=~/.vim/bundle/Vundle.vim
call vundle#begin()

" let Vundle manage Vundle, required
Plugin 'VundleVim/Vundle.vim'
Plugin 'sjl/badwolf'
Plugin 'vim-airline/vim-airline' 
Plugin 'vim-airline/vim-airline-themes'
Plugin 'Valloric/YouCompleteMe'
Plugin 'kien/ctrlp.vim'
Plugin 'scrooloose/nerdtree'
Plugin 'lervag/vimtex'

call vundle#end()            

filetype plugin on
filetype indent on

" " Sets how many lines of history VIM has to remember
set history=700
set autoread
set t_Co=256

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" " => Colors and fonts
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

colorscheme badwolf         " awesome colorscheme
syntax enable           " enable syntax processing

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" " => Tabs and spacing
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
set tabstop=4       " number of visual spaces per TAB
set softtabstop=4   " number of spaces in tab when editing
set expandtab       " tabs are spaces
" Be smart when using tabs ;)
set smarttab
set si "Smart indent
set splitbelow  "Open split windows below
set splitright  "Open split windows to the right
set backspace=indent,eol,start
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" " => UI Config
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
set number              " show line numbers
set showcmd             " show command in bottom bar
"set cursorline          " highlight current line
set wildmenu            " visual autocomplete for command menu
set lazyredraw          " redraw only when we need to.
set showmatch           " highlight matching [{()}]

" Airline
" " Need fugitive to see git branch
set laststatus=2
let g:airline#extensions#tabline#enabled = 1
let g:airline_theme='badwolf'
" the separator used on the left side
let g:airline_left_sep=''
" " the separator used on the right side 
let g:airline_right_sep=''

" Ycm
let g:ycm_extra_conf_globlist = ['~/OpenFOAM/*']
let g:ycm_global_ycm_extra_conf = '~/.vim/ycm_settings/.ycm_extra_conf.py'
let g:ycm_autoclose_preview_window_after_completion = 1
au BufNewFile,BufRead *.C,*.H let g:ycm_filetype_blacklist = {'cpp':1}

" NERDTree
"
nmap <F3> :NERDTreeToggle<CR>

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" " => Visual mode related
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" Visual mode pressing 1 or 2 searches for the current selection
" " Super useful! From an idea by Michael Naumann
vnoremap <silent> 1 :call VisualSelection('f')<CR>
vnoremap <silent> 2 :call VisualSelection('b')<CR>

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" " => Searching
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
set incsearch           " search as characters are entered
set hlsearch            " highlight matches
"" stop highlighting the old search with spacebar

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" " => Folding
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
set foldenable          " enable folding
set foldlevelstart=10   " open most folds by default
set foldnestmax=10      " 10 nested fold max

nnoremap <space> za
set foldmethod=indent   " fold based on indent level

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" " => Movement, tabs, windows and buffers
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

let mapleader=","       " leader is comma

" move vertically by visual line
nnoremap j gj
nnoremap k gk

" move to beginning/end of line
nnoremap b ^
nnoremap e $

set mouse=a ""Enable mouse options

"Remap esc to something closer
inoremap jk <esc>

" highlight last inserted text
nnoremap gv `[v`]

" Smart way to move between windows
map <C-j> <C-W>j
map <C-k> <C-W>k
map <C-h> <C-W>h
map <C-l> <C-W>l

" Close the current buffer
map <leader>bd :Bclose<cr>

" " Close all the buffers
map <leader>ba :1,1000 bd!<cr>

" Useful mappings for managing tabs
map <leader>tt :tabnew<cr>
map <leader>to :tabonly<cr>
map <leader>ty :tabclose<cr>
map <leader>tj :tabnext<cr>
map <leader>tk :tabprev<cr>

nnoremap <leader>v :vsp<Space>
nnoremap <leader>h :split<Space>

"Search for word
vnoremap // y/<C-R>"<CR> 

" Opens a new tab with the current buffer's path
" " Super useful when editing files in the same directory
map <leader>te :tabedit <c-r>=expand("%:p:h")<cr>/

" " Switch CWD to the directory of the open buffer
map <leader>cd :cd %:p:h<cr>:pwd<cr>

" Specify the behavior when switching between buffers 
try
   set switchbuf=useopen,usetab,newtab
   set stal=2
catch
endtry

" Return to last edit position when opening files (You want this!)
autocmd BufReadPost *
     \ if line("'\"") > 0 && line("'\"") <= line("$") |
     \   exe "normal! g`\"" |
     \ endif
" Remember info about open buffers on close
 set viminfo^=%


" Move a line of text using ALT+[jk] or Comamnd+[jk] on mac
nmap <M-j> mz:m+<cr>`z
nmap <M-k> mz:m-2<cr>`z
vmap <M-j> :m'>+<cr>`<my`>mzgv`yo`z
vmap <M-k> :m'<-2<cr>`>my`<mzgv`yo`z

" Delete trailing white space on save, useful for Python
func! DeleteTrailingWS()
   exe "normal mz"
   %s/\s\+$//ge
   exe "normal `z"
endfunc
autocmd BufWrite *.py :call DeleteTrailingWS()

" Quit without saving remapped to q
:cabbrev q q! 

nnoremap <leader>p "0p
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" " => Leader shortcuts
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" toggle gundo: super-undo
nnoremap <leader>u :GundoToggle<CR>

" edit vimrc/zshrc and load vimrc bindings
nnoremap <leader>ev :vsp $MYVIMRC<CR>
nnoremap <leader>ez :vsp ~/.zshrc<CR>
nnoremap <leader>sv :source $MYVIMRC<CR>

" save session
nnoremap <leader>sa :mksession<CR>

""Silver searcher
let g:ackprg = 'ag --vimgrep'
nnoremap <leader>a :Ag

nnoremap <leader>s :mksession!<CR>
nnoremap <leader>g :YcmCompleter GoToImprecise<CR>
nnoremap <leader>f :YcmCompleter FixIt<CR>
nmap     <C-p>     :CtrlP<CR><C-\>w
nnoremap <leader><space> :nohlsearch<CR>
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" " => CtrlP settings
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
set runtimepath^=~/.vim/bundle/ctrlp.vim
let g:ctrlp_match_window = 'bottom,order:ttb'
let g:ctrlp_switch_buffer = 0
let g:ctrlp_working_path_mode = 0
let g:ctrlp_user_command = 'ag %s -l --nocolor --hidden -g ""'


""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" " => Launch config
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"call pathogen#infect()                      " use pathogen
"call pathogen#runtime_append_all_bundles()  " use pathogen


""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" " => Backups
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
set backup
set backupdir=~/.vim-tmp,~/.tmp,~/tmp,/var/tmp,/tmp
set backupskip=/tmp/*,/private/tmp/*
set directory=~/.vim-tmp,~/.tmp,~/tmp,/var/tmp,/tmp
set writebackup

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" " => Spell check
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" Pressing ,ss will toggle and untoggle spell checking
set spelllang=en
map <leader>ss :setlocal spell!<cr>
"
" " Shortcuts using <leader>
map <leader>sn ]s
map <leader>sp [s
map <leader>sa zg
map <leader>sf z=

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" " => Autogroup
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
augroup configgroup
    autocmd!
    autocmd FileType cpp,hpp setlocal tabstop=4
    autocmd FileType cpp,hpp setlocal shiftwidth=4
    autocmd FileType cpp,hpp setlocal softtabstop=4
    autocmd FileType tex setlocal tabstop=2
    autocmd FileType tex setlocal shiftwidth=2
    autocmd FileType tex setlocal softtabstop=2
augroup END

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" " => Custom commands
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"The following command puts comments using # near the left margin 
""when the user types ":Compy" in visual mode.
command Compy :'<,'>s/^/#/
noremap   <buffer> c      :s,^\(\s*\)[^# \t]\@=,\1#,e<CR>:nohls<CR>zvj
noremap   <buffer> <C-c>  :s,^\(\s*\)#\s\@!,\1,e<CR>:nohls<CR>zvj


""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" " => THE END!!
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
