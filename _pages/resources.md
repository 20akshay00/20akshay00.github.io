---
layout: archive
title: "Resources"
permalink: /resources/
author_profile: true
skip_toc: false
---

{% include base_path %}

In today's world, there is no limit to the content that we can access for learning new things. I have often burnt many
days simply collecting and vetting these resources, making it quite hard to make tangible progress in actually learning
the content. This page is dedicated to compiling any useful resources that I have found, in the hopes that someone else
may find it helpful.

<link rel="stylesheet" href="http://yandex.st/highlightjs/6.2/styles/googlecode.min.css">

<script src="http://code.jquery.com/jquery-1.7.2.min.js"></script>
<script src="http://yandex.st/highlightjs/6.2/highlight.min.js"></script>

<script>hljs.initHighlightingOnLoad();</script>
<script type="text/javascript">
    $(document).ready(function () {
        $("h1, h2").each(function (i, item) {
            if ($(this).text() != "Resources")
            {
                var tag = $(item).get(0).localName;
                $(item).attr("id", "section" + i);
                $("#category").append('<a class="new' + tag + '" href="#section' + i + '">' + $(this).text() + '</a></br>');
                $(".newh2").css("margin-left", 20);
            }
        });
    });

</script>

<hr style="text-align:left;margin-left:0;border-top:2px solid #6b7278"> 
<div id="category"></div>
<hr style="text-align:left;margin-left:0;border-top:2px solid #6b7278"> 
# Scientific Programming with Julia
<hr style="text-align:left;margin-left:0;border-top:2px solid #6b7278"> 

{% include_relative /_resources/julia.md %}
<br>

<hr style="text-align:left;margin-left:0;border-top:2px solid #6b7278"> 
# Game Development
<hr style="text-align:left;margin-left:0;border-top:2px solid #6b7278"> 

{% include_relative /_resources/gamedev.md %}
<br>

<hr style="text-align:left;margin-left:0;border-top:2px solid #6b7278"> 
# Visualization
<hr style="text-align:left;margin-left:0;border-top:2px solid #6b7278"> 

{% include_relative /_resources/visualization.md %}
<br>

<hr style="text-align:left;margin-left:0;border-top:2px solid #6b7278"> 
# Quantum Many-Body Physics
<hr style="text-align:left;margin-left:0;border-top:2px solid #6b7278"> 

{% include_relative /_resources/cmt.md %}
<br>
<br>

<hr style="text-align:left;margin-left:0;border-top:2px solid #6b7278"> 
# Random tools
<hr style="text-align:left;margin-left:0;border-top:2px solid #6b7278"> 

{% include_relative /_resources/random_tools.md %}
<hr style="text-align:left;margin-left:0;border-top:2px solid #6b7278"> 

<!-- <ul>
    {% for post in site.resources %}
    <li>{% include archive-single.html %}</li>
    {% endfor %}
</ul> -->