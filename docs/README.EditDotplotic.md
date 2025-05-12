# EditDotplotic.pl

## Overview

`EditDotplotic.pl` is a Perl script that modify an SVG file from Dotplotic

## Features

- modify an SVG file from Dotplotic.

## Usage

```bash
EditDotplotic.pl rotate sbjct_name -30 Dotplotic.svg > Dotplotic.v2.svg
```

### Options

| Option                | Description|Default value|
|-----------------------|------------|-------------|
| `--version`           | Display the program version.| - |
| `--help`              | Display the help message.| - |
| `change`              | Change attributes of the target.| - |
| `delete`              | Delete the target.| - |
| `move`                | Move the target.| - |
| `rotate`              | Rotate the target.| - |


## Requirements

- Perl 5
- Perl modules:
  - `Getopt::Long`
  - `Pod::Text`

## Example

Change the attributes of target:

```bash
EditDotplotic.pl change sbjct_name font-size 7 Dotplotic.svg > Dotplotic.v2.svg
EditDotplotic.pl change sbjct_name text-anchor end Dotplotic.svg > Dotplotic.v2.svg
```

Delete the target:

```bash
EditDotplotic.pl delete sbjct_seq_end > Dotplotic.v2.svg
```

Move the target:

```bash
EditDotplotic.pl move sbjct_name x+30 > Dotplotic.v2.svg
```

Rotate the target:

```bash
EditDotplotic.pl roatate sbjct_name -30 > Dotplotic.v2.svg
```

Multiple actions are separated by ',':

```bash 
EditDotplotic.pl rotate sbjct_name -30, change sbjct_name text-anchor start, change sbjct_name font-size 8, delete *_seq_*, rotate query_name +30, change query_name text-anchor end, change query_name font-size 6 Dotplotic.svg > Dotplotic.v2.svg

```

## License

This script is released under the GNU General Public License v3.0 (GPL-3.0).  
See the [LICENSE](../LICENSE) file for details.
