# Wheat DMC1 and ASY1 ChIP-seq custom scripts

This repository contains custom scripts used in the preparation of the paper entitled "Wheat meiotic recombination hotspots are marked by DMC1 recombinase, the chromosome axis, facultative heterochromatin and signatures of adaptation" by Andrew J. Tock, Daniel M. Holland, Wei Jiang, Kim Osman, Eugenio Sanchez-Moran, James D. Higgins, Keith J. Edwards, Cristobal Uauy, F.C.H. Franklin, and Ian R. Henderson.

These files can be downloaded together by cloning the repository:

```
git clone https://github.com/ajtock/Wheat_DMC1_ASY1_paper
```

## Alignment workflows

Workflows for processing and aligning next-generation sequencing (NGS) reads were developed using [Snakemake](https://snakemake.readthedocs.io/en/stable/) v3.13.3 in conjunction with [conda](https://conda.io/en/latest/) v4.6.9 package and environment manager. These workflows are located in [scripts/read_alignment/](https://github.com/ajtock/Wheat_DMC1_ASY1_paper/tree/master/scripts/read_alignment) in this repository.


Text can be **bold**, _italic_, or ~~strikethrough~~.

## Block quotes

> This is a blockquote following a header.
>
> When something is important enough, you do it even if the odds are not in your favor.

### Code with syntax highlighting

```js
// Javascript code with syntax highlighting.
var fun = function lang(l) {
  dateformat.i18n = require('./lang/' + l)
  return true;
}
```

```ruby
# Ruby code with syntax highlighting
GitHubPages::Dependencies.gems.each do |gem, version|
  s.add_dependency(gem, "= #{version}")
end
```

#### Header 4

*   This is an unordered list following a header.
*   This is an unordered list following a header.
*   This is an unordered list following a header.

##### Header 5

1.  This is an ordered list following a header.
2.  This is an ordered list following a header.
3.  This is an ordered list following a header.

###### Header 6

| head1        | head two          | three |
|:-------------|:------------------|:------|
| ok           | good swedish fish | nice  |
| out of stock | good and plenty   | nice  |
| ok           | good `oreos`      | hmm   |
| ok           | good `zoute` drop | yumm  |

### There's a horizontal rule below this.

* * *

### Here is an unordered list:

*   Item foo
*   Item bar
*   Item baz
*   Item zip

### And an ordered list:

1.  Item one
1.  Item two
1.  Item three
1.  Item four

### And a nested list:

- level 1 item
  - level 2 item
  - level 2 item
    - level 3 item
    - level 3 item
- level 1 item
  - level 2 item
  - level 2 item
  - level 2 item
- level 1 item
  - level 2 item
  - level 2 item
- level 1 item

### Small image

![Octocat](https://github.githubassets.com/images/icons/emoji/octocat.png)

### Large image

![Branching](https://guides.github.com/activities/hello-world/branching.png)


### Definition lists can be used with HTML syntax.

<dl>
<dt>Name</dt>
<dd>Godzilla</dd>
<dt>Born</dt>
<dd>1952</dd>
<dt>Birthplace</dt>
<dd>Japan</dd>
<dt>Color</dt>
<dd>Green</dd>
</dl>

```
Long, single-line code blocks should not wrap. They should horizontally scroll if they are too long. This line should be long enough to demonstrate this.
```

```
The final element.
```
