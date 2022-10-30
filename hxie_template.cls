\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{hxie_template}[2022/10/30 ICPC WF Template By hxie@connect.hku.hk]
\LoadClass[a4paper]{article}

\RequirePackage[margin=5mm, includehead, headsep=0mm, hoffset=0mm]{geometry}

\RequirePackage{fancyhdr}
\pagestyle{fancy}
\rhead{\thepage}

\RequirePackage{color}
\definecolor{dkgreen}{rgb}{0,0.6,0}
\definecolor{gray}{rgb}{0.5,0.5,0.5}
\definecolor{mauve}{rgb}{0.58,0,0.82}

\RequirePackage{titlesec}
\titlespacing{\section}{0pt}{0pt}{0pt}
\titlespacing{\subsection}{0pt}{0pt}{0pt}

\RequirePackage{listings}
\lstset{frame=tb,
  language=c++,
  aboveskip=0mm,
  belowskip=0mm,
  showstringspaces=false,
  columns=flexible,
  basicstyle={\linespread{0.5}\small\ttfamily},
  numbers=none,
  numberstyle=\tiny\color{gray},
  keywordstyle=\color{blue},
  commentstyle=\color{dkgreen},
  stringstyle=\color{mauve},
  breaklines=true,
  breakatwhitespace=true
  tabsize=3
}

\RequirePackage{multicol}
\setlength{\columnseprule}{1pt}

\RequirePackage{graphicx}

\def\uname{No University Name Provided}
\def\tname{No Team Name Provided}
\def\pname{No Participants Name Provided}
\def\lfile{No_Logo_Included}

\newcommand\setinfo{\lhead{\textit{\uname: \textbf{\tname}}}}
\newcommand\university[1]{
    \def\uname{#1}
    \setinfo
}
\newcommand\team[1]{
    \def\tname{#1}
    \setinfo
}
\newcommand\participants[1]{
    \def\pname{#1}
}

\newcommand\ulogo[1]{
    \def\lfile{#1}
}
\newcommand\makecover{
\begin{titlepage}
    \begin{center}
        \vspace*{5cm}
        \Huge
        \textbf{ICPC Team Reference Document} \\
        \vspace*{5cm}
        \tname\\
        \vspace*{1cm}
        \textit{\pname}\\
        \vspace*{5cm}
        \includegraphics[width=3cm]{\lfile}\\
        \vspace*{10pt}
        \textbf{\uname}
    \end{center}
\end{titlepage}
}