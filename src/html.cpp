#include "html.hpp"
using namespace shasta;

#include "iostream.hpp"



void shasta::writeHtmlBegin(ostream& html, const string& title)
{
    html <<
        "<!DOCTYPE html>"
        "<html>"
        "<head>"
        "<meta charset='UTF-8'>"
        "<title>" << title << "</title>";
    writeStyle(html);
    html << "</head>";
}



void shasta::writeHtmlEnd(ostream& html)
{
    html << "</html>";
}



void shasta::writeStyle(ostream& html)
{
    html << R"%(
<style>
    body {
        font-family: Arial;
    }
    pre {
        font-family: courier;
    }
    p, input {
        font-size: 16px;
    }
    h1, h2, h3 {
        color: DarkSlateBlue;
    }
    table {
        border-collapse: collapse;
    }
    th, td {
        border: 1px dashed #b8b5c7d9;
    }
    th {
        font-weight: bold;
        text-align: center;
    }
    th.left {
        text-align: left;
    }
    td.centered {
        text-align: center;
    }
    td.right {
        text-align: right;
    }
    a {
        color: DarkSlateBlue;
    }
    
</style>
    )%";

}
