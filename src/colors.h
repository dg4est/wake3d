/**
 * File:   colors.h
 * Author: akirby
 *
 * Created on October 28, 2024, 2:48 PM
 */

#ifndef COLORS_H
#define COLORS_H

#ifdef __cplusplus
extern "C" {
#endif

/* Reset */
#define COLOR_OFF   "\x1B[0m"         // Text Reset

// Regular Colors
#define BLACK       "\x1B[0;30m"
#define RED         "\x1B[0;31m"
#define GREEN       "\x1B[0;32m"
#define YELLOW      "\x1B[0;33m"
#define BLUE        "\x1B[0;34m"
#define PURPLE      "\x1B[0;35m"
#define CYAN        "\x1B[0;36m"
#define WHITE       "\x1B[0;37m"

// Bold
#define BBLACK      "\x1B[1;30m"
#define BRED        "\x1B[1;31m"
#define BGREEN      "\x1B[1;32m"
#define BYELLOW     "\x1B[1;33m"
#define BBLUE       "\x1B[1;34m"
#define BPURPLE     "\x1B[1;35m"
#define BCYAN       "\x1B[1;36m"
#define BWHITE      "\x1B[1;37m"

// Underline
#define UBLACK      "\x1B[4;30m"
#define URED        "\x1B[4;31m"
#define UGREEN      "\x1B[4;32m"
#define UYELLOW     "\x1B[4;33m"
#define UBLUE       "\x1B[4;34m"
#define UPURPLE     "\x1B[4;35m"
#define UCYAN       "\x1B[4;36m"
#define UWHITE      "\x1B[4;37m"

// Background
#define ON_BLACK    "\x1B[40m"
#define ON_RED      "\x1B[41m"
#define ON_GREEN    "\x1B[42m"
#define ON_YELLOW   "\x1B[43m"
#define ON_BLUE     "\x1B[44m"
#define ON_PURPLE   "\x1B[45m"
#define ON_CYAN     "\x1B[46m"
#define ON_WHITE    "\x1B[47m"

// High Intensity
#define IBLACK      "\x1B[0;90m"
#define IRED        "\x1B[0;91m"
#define IGREEN      "\x1B[0;92m"
#define IYellow     "\x1B[0;93m"
#define IBLUE       "\x1B[0;94m"
#define IPURPLE     "\x1B[0;95m"
#define ICYAN       "\x1B[0;96m"
#define IWHITE      "\x1B[0;97m"

// Bold High Intensity
#define BIBLACK     "\x1B[1;90m"
#define BIRED       "\x1B[1;91m"
#define BIGREEN     "\x1B[1;92m"
#define BIYELLOW    "\x1B[1;93m"
#define BIBLUE      "\x1B[1;94m"
#define BIPURPLE    "\x1B[1;95m"
#define BICYAN      "\x1B[1;96m"
#define BIWHITE     "\x1B[1;97m"

// High Intensity backgrounds
#define ON_IBLACK   "\x1B[0;100m"
#define ON_IRED     "\x1B[0;101m"
#define ON_IGREEN   "\x1B[0;102m"
#define ON_IYELLOW  "\x1B[0;103m"
#define ON_IBLUE    "\x1B[0;104m"
#define ON_IPURPLE  "\x1B[10;95m"
#define ON_ICYAN    "\x1B[0;106m"
#define ON_IWHITE   "\x1B[0;107m"

#ifdef __cplusplus
}
#endif
#endif /* COLORS_H */