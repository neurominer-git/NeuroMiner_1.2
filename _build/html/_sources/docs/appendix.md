(appendix)=
# Appendix


<!-- Figures
=======

This appendix describes how the actual page is build from its components
and how they are influenced by TeX's parameters. (
Figure [\[bild\]](#bild){reference-type="ref" reference="bild"} has been
created by Nelson Beebe at the University of Utah.)

text area

:   The normal text area ("Body") contains the running text including
    footnotes, tables and figures. The headings, footer and margin notes
    do *not* belong to the text area.

    The text area has the width and the height .

    In a two column layout the text area is split into two columns, with
    the width each and a space of between them. Thus the is a little bit
    smaller than half the .

    and should be a multiple of the width of one character in the `tt`
    font.

    should be multiple of the line height , increased by the constant
    value of .

    Indentations inside the text area are defined with and . These
    parameters should not be changed explicitly by the user but rather
    implicitly through environments.

left margin

:   The left margin is either or plus 1 inch. Both parameters have the
    same value, unless the `twoside` option is given.

top margin

:   The top margin is the sum of , and plus 1 inch.

right margin

:   The right margin is the paper width minus the left margin and the
    text area.

bottom margin

:   The bottom margin is the paper height minus the top margin and the
    text area.

heading

:   The heading is inside the top margin with a space of between the
    lower border of the header and the upper border of the text area.
    Above the header is a free space of increased by 1 inch.

footing

:   The footer is inside the bottom margin with a space of between the
    lower border of the text area and the lower border of the footer.

margin notes

:   Margin notes are inside the left or right margin. They have a width
    of and a space of between the margin note and the text area. The
    vertical space between two margin notes is .

The `paperheight` consists of the following elements (from top to
bottom):

> 1 inch\
> \
> \
> \
> \
> \
> remaining page.

On pages with margin notes in the right margin the `paperwidth` consists
of the following elements:

> 1 inch\
>   or  \
> \
> \
> \
> remaining page

With the option `twoside` the left pages change to

> 1 inch\
> \
> \
> remaining page

The parameters , , and may be negative. In this case, the margin will be
smaller than 1 inch. The same is true for and which leads to text that
is wider than the text area.

Extensive treatment and figures to this subject may be found in the
[TUGboat]{.smallcaps} Vol.9, No.1 (April 1988).

The parameter is no longer defined in  since no-one used it.

by -1000 by 2

by 2 by 3613

= by 3613 by 1000

by 50400 by 1807 by 2400 by 1807 by 50400 by 4207

by 32522 by 10841

by 32522 by 2800

by 75 by 100 by 4207

= by -4444 by -2000

= by 4444

by -7227

= by -

by 32522

by 4 by 10

by 3

(61430,79497)

(0,0)

(0,0)(-0,-0) (61430,79497)\[\]

(0,)

(0,0)(-0,-0) (0,0)(1,0)7227 (7227,0)(-1,0)7227

(7227,)

(0,0)(-0,-0) (0,0)(0,1)7227 (0,7227)(0,-1)7227

(7227,)

(,) (0,0)(0,1) (0,)(1,0) (,0)(0,1) (3613,)

(0,0)(-0,-0) (32522, 2000)\[l\]A line of text...

(3613,)

(0,0)(-0,--2000) ( 32522,2000)\[l\]Next line...

(,)

(0,0)(-0,--2000) (0,0)(0,1) 2000 (0, 2000)(0,-1) 2000

(0,50400)

(0,0)(-0,-0) (0,0)(1,0)3613 (3613,0)(-1,0)3613

(0,50400)

(0,0)(-0,-0) (0,0)(1,0)3613 (3613,0)(-1,0)3613

(3613,4207)

(0,0)(-0,-0) (32522, 50400)\[\]Page Text

(3613,0)

(0,0)(-0,-0) (32522, 2400)\[\]Page Footer

(3613,0)

(0,0)(--800,-0) (0,0)(0,1)2400 (0,2400)(0,-1)2400

(,0)

(0,0)(-800,-0) (0,0)(0,1)4207 (0,4207)(0,-1)4207

(3613,)

(0,0)(-0,-0) (0,0)(1,0)32522 (32522,0)(-1,0)32522

(3613,4207)

(0,0)(--800,-0) (0,0)(0,1) 50400 (0, 50400)(0,-1) 50400

(3613,)

(0,0)(--800,-0) (0,0)(0,1) 2400 (0, 2400)(0,-1) 2400

(,)

(0,0)(-800,--1807) (0,0)(0,1) 1807 (0, 1807)(0,-1) 1807

(,)

(0,0)(-800,--5420) (0,0)(0,1) 5420 (0, 5420)(0,-1) 5420

(3613,)

(0,0)(-0,-0) (32522, 2400)\[\]Page Header

(,)

(0,0)(-0,-0) ( 7227,4444)\[\]

(,)

(0,0)(--800,- -2000) (0,0)(0,1)2000 (0,2000)(0,-1)2000

(,)

(0,0)(-0,-0) ( 7227,4444)\[\]

(,)

(0,0)(-0,--800) (0,0)(1,0) 7227 ( 7227,0)(-1,0) 7227

(,)

(0,0)(- -2800,-800) (0,0)(1,0)2800 (2800,0)(-1,0)2800

Tables
======

s[\[refman\]]{#refman label="refman"}

The `refman.sty` was defined at the EDV-Zentrum (computing center) of
the TU[^1] Wien. This layout is suitable for reference manuals,
technical descriptions and similar applications. It is based on the
ideas shown in previous sections: The layout has a wide left margin for
headings and margin notes and smaller margins on the right side, the top
and the bottom.

In 1994 this layout was re-implemented as a class for the new . This
made it possible to include some minor improvements, such as the support
of different paper sizes. The `refman.sty` was split into two classes
`refrep`, similar to `report` and `refart`, similar to `article`. These
classes differ in the layout of the header and footer. The `refart` does
not support the command.

The current version of both classes is described in this document. It
serves as an example for the layout.

Invocation
----------

The LaTeX local guide (if available) shows if this class is available at
your TeX installation or where to install it. To use the `refart` class,
simply call it with the command:

    \documentclass[11pt,a4paper]{refart}
    \usepackage{german} % other packages you may want

Options
-------

The `refart` class replaces `article` and `refrep` replaces `report`.
They support all options of these classes except for the `twocolumn`
option.

It supports the additional option `square` which makes the equal to the
.

Neither `refart` nor `refrep` support two column layout, thus the
commands and must not be used.

The index will be set in two column format and you can't change it with
the means of this class.

Layout changes
--------------

### Page design

In this design the usable area for text () is calculated as the paper
width minus 2 . The default value for is 1 Inch.

The option `smallborder` reduces to 0.25 Inch. This is more suitable for
documents viewed on screen, especially when combined with the `a5paper`
and `landscape` options.

Only a fraction of this width is used for the running text (), the
remaining part forms a wide left margin () which is used for headings
and margin notes. The is 70 % of the by default, but this can be changed
with the command which accepts arguments between 0 and 1.

The text height is calculated as the paper height minus 2 . The
`topmargin` is modified by some pagestyles. (see
[2.4.3](#pagestyle){reference-type="ref" reference="pagestyle"}).

The pages are always set with a ragged bottom.

### Section headings

The headings for , , and extend into the left margin, thus using the
full width of the page. They are not justified and hyphenation is
discouraged. A small space is kept free above and below the heading.
Headings for and are set in a bold font.

The `refrep` class defines a different layout for the command: It always
starts a new page and prints the chapter headings in a large bold font
with a thick line above and below. This heading uses the full width of
the page.

A similar heading is created by the commands which is available in both
classes. It uses a roman part number instead of the arabic section
number.

The commands sets the title of the document in the same layout when no
special title page is requested. (This is the default for `refart`. To
suppress the title page in a `refrep` document, you can use the
`notitlepage` option.) The name of the author and the date is printed in
italic flush right below the document title.

### Paragraphs

Paragraphs are separated by a vertical space () of half a line
(`0.5\baselineskip`) plus a stretchable length of 2 pt. Paragraphs are
not indented.

The vertical spacing inside, above and below a list environment is the
same as in the running text.

Footnotes
---------

The footnote layout consists of a small margin (1em) which contains the
footnote symbol. A small space is set between the symbol and the
footnote text. The paragraphs of the footnote are not indented. There is
currently no space between two footnotes, I'm not sure it this will stay
this way. The footnote symbol is set as a superscript. This may change
in later versions. I'm relying on user feedback to finally solve this.

### Description environment

The `description` environment will use the whole left margin for the
description label.

You will find examples in the
section [\[layout\]](#layout){reference-type="ref" reference="layout"}.

### Positioning of margin notes

Margin notes () are always put into the left margin. They use the whole
width of the margin.

The minimum space between two margin notes is set to 0 to prevent them
from being shifted around when many margin notes are used.

### Headers and Footers {#pagestyle}

The page style `plain` puts the page number into the footer in the right
corner. When the option `twoside` is active, the page number of left
pages is put into the left corner.

The pagestyles `headings` and `myheadings` create a header which spans
the whole width of the page. The headings contain the running head ( and
in `refart` and and in `refrep`) when `headings` is used or a fixed text
that can defined with the command when `myheadings` is used. The heading
will be set in a slanted font and separated from the body by a thin
line.

In addition to the standard classes, `refman` supports a style for
footers, which is used in this documentation. The information is exactly
the same as in the headings but now printed in the footer with a thin
line above.

To use a user-defined string you can say:

    \pagestyle{myfootings} % or myheadings
    \markboth{left title}{right title}

The `heading` and `myheading` commands increase the top margin by one
line while the `footings` and `myfootings` commands decrease the top
margin by one line. The page styles `empty` and `plain` leave the top
margin unchanged. You should not combine headings and footings in one
document.

User feedback has shown that it is not a good idea to combine `plain`
and `(my)heading` either. Therefore I changed the layout of the page to
`empty`. Maybe it is necessary to define a `hplain` and `fplain`
pagestyle or to define some magic to use the correct definition of
`plain`. Feedback is welcome.

Additional commands
-------------------

### Marginlabel

The command `\marginlabel{xxx}` prints the text `xxx` right justified
into the left margin. Please note that a will print it left justified.

The word "Example" in the left margin is printed with the command
`\marginlabel{Example:}`

### Attention

The command puts an exclamation mark with an arrow pointing to the text
into the left margin. This is an example for .

Since version 2.0c you can change the symbol used for the command using
a `\renewcommand{\attentionsymbol}{\texttt\{:-)\}}` command. To get the
default back use `\renewcommand{\attentionsymbol}`
`{\large \bfseries ! $\rightarrow$}`

Since version 2.0c takes an optional argument to define the symbol used
in the margin. Thus you can change the symbol once, without having to
restore it later. Do not abuse this feature, it is primarily meant as an
support for the `manfnt` package which enables you to use the "dangerous
bend" and "double dangerous bend" signs.

The `manfnt` package is no longer enclosed with Refman, it has grown and
is now a package of its own.

### Seealso

The command `\seealso{n}` prints an arrow and its argument into the left
margin. You will find examples for this in the left margin and in
chapter 1.

### Maxipage environment

The `maxipage` environment is a special kind of `minipage` which extends
over the full width of the page. It can be used for long formulas or
`tabular` environments. You may use `maxipage` environments inside
floats. You cannot use margin notes inside a `maxipage` and no page
break will occur while in a `maxipage`. A `maxipage` is always a
paragraph of its own with a thick line above and below. You can disable
these lines with the command. They are on by default.

The following paragraph is an example for a `maxipage`:

This very long line is an example for a `maxipage`. It extends over the
full width of the page, including the left margin.

This is normal text after the `maxipage`.

### Fullpage environment {#fullpage}

The `fullpage` environment consists of one or more pages where the text
extends over the full width of the page. You cannot use margin notes
inside a `fullpage` environment. A `fullpage` will always start and end
on a page of its own. It may be used for large tables, program listings
or anything that does not fit into the normal layout.

Page  is an example for a `fullpage`.

### Noparskip

The removes the vertical space between two paragraphs. It is similar to
the command that removes the indent of the first line of a paragraph.

### Setleftmarginwidth

The command is no longer supported. You can achieve similar results by
using the command.

### Descriptioncolon

By default a colon is printed after the description label. The command
disables the colon, the re-enables it.

### Descriptionleft

The command sets the description label left justified into the margin.
The default is right justified which will be achieved with

### Maxipagerule

You can disable the rules before and after a `maxipage` with the command
and re-enable them with the command. The default is on. You should not
mix `maxipage`s with and without rules in one document.

### Condbreak

The command `\condbreak{2cm}` ensures, that the next 2 cm are either
completely on this page or completely on the next. No page break will
appear in the next 2 cm.

This is really a hack to achieve what the command often fails to do.

### Example

The `example` environment acts like a `verse` environment but uses a
`tt` font.

### Pageperchapter

The command creates page number that start with 1 for every new chapter.
This may be useful for larger manuals. Since it works with chapters it
is only available in the `refrep` class.

### Smallborder

The normal border around the page is 1 Inch. That is fine for a printed
document, but wastes a lot of space when a document is meant for reading
on screen. The option `smallborder` reduces the margin to 0.25 Inch.

You can redefine the border with `{0.25in}`. Call
`\setpagefraction{0.7}` afterwards to recalculate the page layout.

### Dvips

The option `dvips` tells DVIPS about the current page size.

### Pdftex

The option `pdftex` tells PDFTeXabout the current page size.

### Pagesize

`pagesize` chooses the correct -command to tell the DVI-driver about the
paper size. It works with DVIPS and DVIPDFMX for DVI output and
PDFTeX for PDF output.

### Ifpdfoutput

You can use `\ifpdfoutput{pdftext}{dvitext}` to write different text
depending on the output format. This command was necessary to implement
the `pagesize` option and is available for the user as well.

The last four commands have been taken from KOMA-Script, thanks Markus.

[^1]: Technical University -->
