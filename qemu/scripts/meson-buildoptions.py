#! /usr/bin/env python3

# Generate configure command line options handling code, based on Meson's
# user build options introspection data
#
# Copyright (C) 2021 Red Hat, Inc.
#
# Author: Paolo Bonzini <pbonzini@redhat.com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import json
import textwrap
import shlex
import sys

SKIP_OPTIONS = {
    "audio_drv_list",
    "default_devices",
    "docdir",
    "fuzzing_engine",
    "qemu_firmwarepath",
    "qemu_suffix",
    "sphinx_build",
    "trace_file",
}

LINE_WIDTH = 76


# Convert the default value of an option to the string used in
# the help message
def value_to_help(value):
    if isinstance(value, list):
        return ",".join(value)
    if isinstance(value, bool):
        return "enabled" if value else "disabled"
    return str(value)


def wrap(left, text, indent):
    spaces = " " * indent
    if len(left) >= indent:
        yield left
        left = spaces
    else:
        left = (left + spaces)[0:indent]
    yield from textwrap.wrap(
        text, width=LINE_WIDTH, initial_indent=left, subsequent_indent=spaces
    )


def sh_print(line=""):
    print('  printf "%s\\n"', shlex.quote(line))


def help_line(left, opt, indent, long):
    right = f'{opt["description"]}'
    if long:
        value = value_to_help(opt["value"])
        if value != "auto":
            right += f" [{value}]"
    if "choices" in opt and long:
        choices = "/".join(sorted(opt["choices"]))
        right += f" (choices: {choices})"
    for x in wrap("  " + left, right, indent):
        sh_print(x)


# Return whether the option (a dictionary) can be used with
# arguments.  Booleans can never be used with arguments;
# combos allow an argument only if they accept other values
# than "auto", "enabled", and "disabled".
def allow_arg(opt):
    if opt["type"] == "boolean":
        return False
    if opt["type"] != "combo":
        return True
    return not (set(opt["choices"]) <= {"auto", "disabled", "enabled"})


def load_options(json):
    json = [
        x
        for x in json
        if x["section"] == "user"
        and ":" not in x["name"]
        and x["name"] not in SKIP_OPTIONS
    ]
    return sorted(json, key=lambda x: x["name"])


def print_help(options):
    print("meson_options_help() {")
    for opt in options:
        key = opt["name"].replace("_", "-")
        # The first section includes options that have an arguments,
        # and booleans (i.e., only one of enable/disable makes sense)
        if opt["type"] == "boolean":
            left = f"--disable-{key}" if opt["value"] else f"--enable-{key}"
            help_line(left, opt, 27, False)
        elif allow_arg(opt):
            if opt["type"] == "combo" and "enabled" in opt["choices"]:
                left = f"--enable-{key}[=CHOICE]"
            else:
                left = f"--enable-{key}=CHOICE"
            help_line(left, opt, 27, True)

    sh_print()
    sh_print("Optional features, enabled with --enable-FEATURE and")
    sh_print("disabled with --disable-FEATURE, default is enabled if available")
    sh_print("(unless built with --without-default-features):")
    sh_print()
    for opt in options:
        key = opt["name"].replace("_", "-")
        if opt["type"] != "boolean" and not allow_arg(opt):
            help_line(key, opt, 18, False)
    print("}")


def print_parse(options):
    print("_meson_option_parse() {")
    print("  case $1 in")
    for opt in options:
        key = opt["name"].replace("_", "-")
        name = opt["name"]
        if opt["type"] == "boolean":
            print(f'    --enable-{key}) printf "%s" -D{name}=true ;;')
            print(f'    --disable-{key}) printf "%s" -D{name}=false ;;')
        else:
            if opt["type"] == "combo" and "enabled" in opt["choices"]:
                print(f'    --enable-{key}) printf "%s" -D{name}=enabled ;;')
            if opt["type"] == "combo" and "disabled" in opt["choices"]:
                print(f'    --disable-{key}) printf "%s" -D{name}=disabled ;;')
            if allow_arg(opt):
                print(f'    --enable-{key}=*) quote_sh "-D{name}=$2" ;;')
    print("    *) return 1 ;;")
    print("  esac")
    print("}")


def fixup_options(options):
    # Meson <= 0.60 does not include the choices in array options, fix that up
    for opt in options:
        if opt["name"] == "trace_backends":
            opt["choices"] = [
                "dtrace",
                "ftrace",
                "log",
                "nop",
                "simple",
                "syslog",
                "ust",
            ]


options = load_options(json.load(sys.stdin))
fixup_options(options)
print("# This file is generated by meson-buildoptions.py, do not edit!")
print_help(options)
print_parse(options)
