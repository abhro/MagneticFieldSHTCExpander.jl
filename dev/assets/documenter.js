// Generated by Documenter.jl
requirejs.config({
  paths: {
    'highlight-julia': 'https://cdnjs.cloudflare.com/ajax/libs/highlight.js/11.8.0/languages/julia.min',
    'headroom': 'https://cdnjs.cloudflare.com/ajax/libs/headroom/0.12.0/headroom.min',
    'jqueryui': 'https://cdnjs.cloudflare.com/ajax/libs/jqueryui/1.13.2/jquery-ui.min',
    'katex-auto-render': 'https://cdnjs.cloudflare.com/ajax/libs/KaTeX/0.16.8/contrib/auto-render.min',
    'jquery': 'https://cdnjs.cloudflare.com/ajax/libs/jquery/3.7.0/jquery.min',
    'headroom-jquery': 'https://cdnjs.cloudflare.com/ajax/libs/headroom/0.12.0/jQuery.headroom.min',
    'katex': 'https://cdnjs.cloudflare.com/ajax/libs/KaTeX/0.16.8/katex.min',
    'highlight': 'https://cdnjs.cloudflare.com/ajax/libs/highlight.js/11.8.0/highlight.min',
    'highlight-julia-repl': 'https://cdnjs.cloudflare.com/ajax/libs/highlight.js/11.8.0/languages/julia-repl.min',
  },
  shim: {
  "highlight-julia": {
    "deps": [
      "highlight"
    ]
  },
  "katex-auto-render": {
    "deps": [
      "katex"
    ]
  },
  "headroom-jquery": {
    "deps": [
      "jquery",
      "headroom"
    ]
  },
  "highlight-julia-repl": {
    "deps": [
      "highlight"
    ]
  }
}
});
////////////////////////////////////////////////////////////////////////////////
require(['jquery', 'katex', 'katex-auto-render'], function($, katex, renderMathInElement) {
$(document).ready(function() {
  renderMathInElement(
    document.body,
    {
  "delimiters": [
    {
      "left": "$",
      "right": "$",
      "display": false
    },
    {
      "left": "$$",
      "right": "$$",
      "display": true
    },
    {
      "left": "\\[",
      "right": "\\]",
      "display": true
    }
  ],
  "macros": {
    "\\grad": "\\boldsymbol{∇}",
    "\\dpdv": "\\dfrac{∂#1}{∂#2}",
    "\\pdv": "\\frac{∂#1}{∂#2}",
    "\\TT": "\\mathsf{T}",
    "\\dv": "\\frac{d#1}{d#2}",
    "\\ddv": "\\dfrac{d#1}{d#2}"
  }
}

  );
})

})
////////////////////////////////////////////////////////////////////////////////
require(['jquery', 'highlight', 'highlight-julia', 'highlight-julia-repl'], function($) {
$(document).ready(function() {
    hljs.highlightAll();
})

})
////////////////////////////////////////////////////////////////////////////////
require(['jquery'], function($) {

let timer = 0;
var isExpanded = true;

$(document).on(
  "click",
  ".docstring .docstring-article-toggle-button",
  function () {
    let articleToggleTitle = "Expand docstring";
    const parent = $(this).parent();

    debounce(() => {
      if (parent.siblings("section").is(":visible")) {
        parent
          .find("a.docstring-article-toggle-button")
          .removeClass("fa-chevron-down")
          .addClass("fa-chevron-right");
      } else {
        parent
          .find("a.docstring-article-toggle-button")
          .removeClass("fa-chevron-right")
          .addClass("fa-chevron-down");

        articleToggleTitle = "Collapse docstring";
      }

      parent
        .children(".docstring-article-toggle-button")
        .prop("title", articleToggleTitle);
      parent.siblings("section").slideToggle();
    });
  }
);

$(document).on("click", ".docs-article-toggle-button", function (event) {
  let articleToggleTitle = "Expand docstring";
  let navArticleToggleTitle = "Expand all docstrings";
  let animationSpeed = event.noToggleAnimation ? 0 : 400;

  debounce(() => {
    if (isExpanded) {
      $(this).removeClass("fa-chevron-up").addClass("fa-chevron-down");
      $("a.docstring-article-toggle-button")
        .removeClass("fa-chevron-down")
        .addClass("fa-chevron-right");

      isExpanded = false;

      $(".docstring section").slideUp(animationSpeed);
    } else {
      $(this).removeClass("fa-chevron-down").addClass("fa-chevron-up");
      $("a.docstring-article-toggle-button")
        .removeClass("fa-chevron-right")
        .addClass("fa-chevron-down");

      isExpanded = true;
      articleToggleTitle = "Collapse docstring";
      navArticleToggleTitle = "Collapse all docstrings";

      $(".docstring section").slideDown(animationSpeed);
    }

    $(this).prop("title", navArticleToggleTitle);
    $(".docstring-article-toggle-button").prop("title", articleToggleTitle);
  });
});

function debounce(callback, timeout = 300) {
  if (Date.now() - timer > timeout) {
    callback();
  }

  clearTimeout(timer);

  timer = Date.now();
}

})
////////////////////////////////////////////////////////////////////////////////
require([], function() {
function addCopyButtonCallbacks() {
  for (const el of document.getElementsByTagName("pre")) {
    const button = document.createElement("button");
    button.classList.add("copy-button", "fa-solid", "fa-copy");
    button.setAttribute("aria-label", "Copy this code block");
    button.setAttribute("title", "Copy");

    el.appendChild(button);

    const success = function () {
      button.classList.add("success", "fa-check");
      button.classList.remove("fa-copy");
    };

    const failure = function () {
      button.classList.add("error", "fa-xmark");
      button.classList.remove("fa-copy");
    };

    button.addEventListener("click", function () {
      copyToClipboard(el.innerText).then(success, failure);

      setTimeout(function () {
        button.classList.add("fa-copy");
        button.classList.remove("success", "fa-check", "fa-xmark");
      }, 5000);
    });
  }
}

function copyToClipboard(text) {
  // clipboard API is only available in secure contexts
  if (window.navigator && window.navigator.clipboard) {
    return window.navigator.clipboard.writeText(text);
  } else {
    return new Promise(function (resolve, reject) {
      try {
        const el = document.createElement("textarea");
        el.textContent = text;
        el.style.position = "fixed";
        el.style.opacity = 0;
        document.body.appendChild(el);
        el.select();
        document.execCommand("copy");

        resolve();
      } catch (err) {
        reject(err);
      } finally {
        document.body.removeChild(el);
      }
    });
  }
}

if (document.readyState === "loading") {
  document.addEventListener("DOMContentLoaded", addCopyButtonCallbacks);
} else {
  addCopyButtonCallbacks();
}

})
////////////////////////////////////////////////////////////////////////////////
require(['jquery', 'headroom', 'headroom-jquery'], function($, Headroom) {

// Manages the top navigation bar (hides it when the user starts scrolling down on the
// mobile).
window.Headroom = Headroom; // work around buggy module loading?
$(document).ready(function () {
  $("#documenter .docs-navbar").headroom({
    tolerance: { up: 10, down: 10 },
  });
});

})
////////////////////////////////////////////////////////////////////////////////
require(['jquery'], function($) {

$(document).ready(function () {
  let meta = $("div[data-docstringscollapsed]").data();

  if (meta?.docstringscollapsed) {
    $("#documenter-article-toggle-button").trigger({
      type: "click",
      noToggleAnimation: true,
    });
  }
});

})
////////////////////////////////////////////////////////////////////////////////
require(['jquery'], function($) {

/*
To get an in-depth about the thought process you can refer: https://hetarth02.hashnode.dev/series/gsoc

PSEUDOCODE:

Searching happens automatically as the user types or adjusts the selected filters.
To preserve responsiveness, as much as possible of the slow parts of the search are done
in a web worker. Searching and result generation are done in the worker, and filtering and
DOM updates are done in the main thread. The filters are in the main thread as they should
be very quick to apply. This lets filters be changed without re-searching with minisearch
(which is possible even if filtering is on the worker thread) and also lets filters be
changed _while_ the worker is searching and without message passing (neither of which are
possible if filtering is on the worker thread)

SEARCH WORKER:

Import minisearch

Build index

On message from main thread
  run search
  find the first 200 unique results from each category, and compute their divs for display
    note that this is necessary and sufficient information for the main thread to find the
    first 200 unique results from any given filter set
  post results to main thread

MAIN:

Launch worker

Declare nonconstant globals (worker_is_running,  last_search_text, unfiltered_results)

On text update
  if worker is not running, launch_search()

launch_search
  set worker_is_running to true, set last_search_text to the search text
  post the search query to worker

on message from worker
  if last_search_text is not the same as the text in the search field,
    the latest search result is not reflective of the latest search query, so update again
    launch_search()
  otherwise
    set worker_is_running to false

  regardless, display the new search results to the user
  save the unfiltered_results as a global
  update_search()

on filter click
  adjust the filter selection
  update_search()

update_search
  apply search filters by looping through the unfiltered_results and finding the first 200
    unique results that match the filters

  Update the DOM
*/

/////// SEARCH WORKER ///////

function worker_function(documenterSearchIndex, documenterBaseURL, filters) {
  importScripts(
    "https://cdn.jsdelivr.net/npm/minisearch@6.1.0/dist/umd/index.min.js"
  );

  let data = documenterSearchIndex.map((x, key) => {
    x["id"] = key; // minisearch requires a unique for each object
    return x;
  });

  // list below is the lunr 2.1.3 list minus the intersect with names(Base)
  // (all, any, get, in, is, only, which) and (do, else, for, let, where, while, with)
  // ideally we'd just filter the original list but it's not available as a variable
  const stopWords = new Set([
    "a",
    "able",
    "about",
    "across",
    "after",
    "almost",
    "also",
    "am",
    "among",
    "an",
    "and",
    "are",
    "as",
    "at",
    "be",
    "because",
    "been",
    "but",
    "by",
    "can",
    "cannot",
    "could",
    "dear",
    "did",
    "does",
    "either",
    "ever",
    "every",
    "from",
    "got",
    "had",
    "has",
    "have",
    "he",
    "her",
    "hers",
    "him",
    "his",
    "how",
    "however",
    "i",
    "if",
    "into",
    "it",
    "its",
    "just",
    "least",
    "like",
    "likely",
    "may",
    "me",
    "might",
    "most",
    "must",
    "my",
    "neither",
    "no",
    "nor",
    "not",
    "of",
    "off",
    "often",
    "on",
    "or",
    "other",
    "our",
    "own",
    "rather",
    "said",
    "say",
    "says",
    "she",
    "should",
    "since",
    "so",
    "some",
    "than",
    "that",
    "the",
    "their",
    "them",
    "then",
    "there",
    "these",
    "they",
    "this",
    "tis",
    "to",
    "too",
    "twas",
    "us",
    "wants",
    "was",
    "we",
    "were",
    "what",
    "when",
    "who",
    "whom",
    "why",
    "will",
    "would",
    "yet",
    "you",
    "your",
  ]);

  let index = new MiniSearch({
    fields: ["title", "text"], // fields to index for full-text search
    storeFields: ["location", "title", "text", "category", "page"], // fields to return with results
    processTerm: (term) => {
      let word = stopWords.has(term) ? null : term;
      if (word) {
        // custom trimmer that doesn't strip @ and !, which are used in julia macro and function names
        word = word
          .replace(/^[^a-zA-Z0-9@!]+/, "")
          .replace(/[^a-zA-Z0-9@!]+$/, "");

        word = word.toLowerCase();
      }

      return word ?? null;
    },
    // add . as a separator, because otherwise "title": "Documenter.Anchors.add!", would not
    // find anything if searching for "add!", only for the entire qualification
    tokenize: (string) => string.split(/[\s\-\.]+/),
    // options which will be applied during the search
    searchOptions: {
      prefix: true,
      boost: { title: 100 },
      fuzzy: 2,
    },
  });

  index.addAll(data);

  /**
   *  Used to map characters to HTML entities.
   * Refer: https://github.com/lodash/lodash/blob/main/src/escape.ts
   */
  const htmlEscapes = {
    "&": "&amp;",
    "<": "&lt;",
    ">": "&gt;",
    '"': "&quot;",
    "'": "&#39;",
  };

  /**
   * Used to match HTML entities and HTML characters.
   * Refer: https://github.com/lodash/lodash/blob/main/src/escape.ts
   */
  const reUnescapedHtml = /[&<>"']/g;
  const reHasUnescapedHtml = RegExp(reUnescapedHtml.source);

  /**
   * Escape function from lodash
   * Refer: https://github.com/lodash/lodash/blob/main/src/escape.ts
   */
  function escape(string) {
    return string && reHasUnescapedHtml.test(string)
      ? string.replace(reUnescapedHtml, (chr) => htmlEscapes[chr])
      : string || "";
  }

  /**
   * RegX escape function from MDN
   * Refer: https://developer.mozilla.org/en-US/docs/Web/JavaScript/Guide/Regular_Expressions#escaping
   */
  function escapeRegExp(string) {
    return string.replace(/[.*+?^${}()|[\]\\]/g, "\\$&"); // $& means the whole matched string
  }

  /**
   * Make the result component given a minisearch result data object and the value
   * of the search input as queryString. To view the result object structure, refer:
   * https://lucaong.github.io/minisearch/modules/_minisearch_.html#searchresult
   *
   * @param {object} result
   * @param {string} querystring
   * @returns string
   */
  function make_search_result(result, querystring) {
    let search_divider = `<div class="search-divider w-100"></div>`;
    let display_link =
      result.location.slice(Math.max(0), Math.min(50, result.location.length)) +
      (result.location.length > 30 ? "..." : ""); // To cut-off the link because it messes with the overflow of the whole div

    if (result.page !== "") {
      display_link += ` (${result.page})`;
    }
    searchstring = escapeRegExp(querystring);
    let textindex = new RegExp(`${searchstring}`, "i").exec(result.text);
    let text =
      textindex !== null
        ? result.text.slice(
            Math.max(textindex.index - 100, 0),
            Math.min(
              textindex.index + querystring.length + 100,
              result.text.length
            )
          )
        : ""; // cut-off text before and after from the match

    text = text.length ? escape(text) : "";

    let display_result = text.length
      ? "..." +
        text.replace(
          new RegExp(`${escape(searchstring)}`, "i"), // For first occurrence
          '<span class="search-result-highlight py-1">$&</span>'
        ) +
        "..."
      : ""; // highlights the match

    let in_code = false;
    if (!["page", "section"].includes(result.category.toLowerCase())) {
      in_code = true;
    }

    // We encode the full url to escape some special characters which can lead to broken links
    let result_div = `
        <a href="${encodeURI(
          documenterBaseURL + "/" + result.location
        )}" class="search-result-link w-100 is-flex is-flex-direction-column gap-2 px-4 py-2">
          <div class="w-100 is-flex is-flex-wrap-wrap is-justify-content-space-between is-align-items-flex-start">
            <div class="search-result-title has-text-weight-bold ${
              in_code ? "search-result-code-title" : ""
            }">${escape(result.title)}</div>
            <div class="property-search-result-badge">${result.category}</div>
          </div>
          <p>
            ${display_result}
          </p>
          <div
            class="has-text-left"
            style="font-size: smaller;"
            title="${result.location}"
          >
            <i class="fas fa-link"></i> ${display_link}
          </div>
        </a>
        ${search_divider}
      `;

    return result_div;
  }

  self.onmessage = function (e) {
    let query = e.data;
    let results = index.search(query, {
      filter: (result) => {
        // Only return relevant results
        return result.score >= 1;
      },
      combineWith: "AND",
    });

    // Pre-filter to deduplicate and limit to 200 per category to the extent
    // possible without knowing what the filters are.
    let filtered_results = [];
    let counts = {};
    for (let filter of filters) {
      counts[filter] = 0;
    }
    let present = {};

    for (let result of results) {
      cat = result.category;
      cnt = counts[cat];
      if (cnt < 200) {
        id = cat + "---" + result.location;
        if (present[id]) {
          continue;
        }
        present[id] = true;
        filtered_results.push({
          location: result.location,
          category: cat,
          div: make_search_result(result, query),
        });
      }
    }

    postMessage(filtered_results);
  };
}

/////// SEARCH MAIN ///////

function runSearchMainCode() {
  // `worker = Threads.@spawn worker_function(documenterSearchIndex)`, but in JavaScript!
  const filters = [
    ...new Set(documenterSearchIndex["docs"].map((x) => x.category)),
  ];
  const worker_str =
    "(" +
    worker_function.toString() +
    ")(" +
    JSON.stringify(documenterSearchIndex["docs"]) +
    "," +
    JSON.stringify(documenterBaseURL) +
    "," +
    JSON.stringify(filters) +
    ")";
  const worker_blob = new Blob([worker_str], { type: "text/javascript" });
  const worker = new Worker(URL.createObjectURL(worker_blob));

  // Whether the worker is currently handling a search. This is a boolean
  // as the worker only ever handles 1 or 0 searches at a time.
  var worker_is_running = false;

  // The last search text that was sent to the worker. This is used to determine
  // if the worker should be launched again when it reports back results.
  var last_search_text = "";

  // The results of the last search. This, in combination with the state of the filters
  // in the DOM, is used compute the results to display on calls to update_search.
  var unfiltered_results = [];

  // Which filter is currently selected
  var selected_filter = "";

  $(document).on("input", ".documenter-search-input", function (event) {
    if (!worker_is_running) {
      launch_search();
    }
  });

  function launch_search() {
    worker_is_running = true;
    last_search_text = $(".documenter-search-input").val();
    worker.postMessage(last_search_text);
  }

  worker.onmessage = function (e) {
    if (last_search_text !== $(".documenter-search-input").val()) {
      launch_search();
    } else {
      worker_is_running = false;
    }

    unfiltered_results = e.data;
    update_search();
  };

  $(document).on("click", ".search-filter", function () {
    if ($(this).hasClass("search-filter-selected")) {
      selected_filter = "";
    } else {
      selected_filter = $(this).text().toLowerCase();
    }

    // This updates search results and toggles classes for UI:
    update_search();
  });

  /**
   * Make/Update the search component
   */
  function update_search() {
    let querystring = $(".documenter-search-input").val();

    if (querystring.trim()) {
      if (selected_filter == "") {
        results = unfiltered_results;
      } else {
        results = unfiltered_results.filter((result) => {
          return selected_filter == result.category.toLowerCase();
        });
      }

      let search_result_container = ``;
      let modal_filters = make_modal_body_filters();
      let search_divider = `<div class="search-divider w-100"></div>`;

      if (results.length) {
        let links = [];
        let count = 0;
        let search_results = "";

        for (var i = 0, n = results.length; i < n && count < 200; ++i) {
          let result = results[i];
          if (result.location && !links.includes(result.location)) {
            search_results += result.div;
            count++;
            links.push(result.location);
          }
        }

        if (count == 1) {
          count_str = "1 result";
        } else if (count == 200) {
          count_str = "200+ results";
        } else {
          count_str = count + " results";
        }
        let result_count = `<div class="is-size-6">${count_str}</div>`;

        search_result_container = `
              <div class="is-flex is-flex-direction-column gap-2 is-align-items-flex-start">
                  ${modal_filters}
                  ${search_divider}
                  ${result_count}
                  <div class="is-clipped w-100 is-flex is-flex-direction-column gap-2 is-align-items-flex-start has-text-justified mt-1">
                    ${search_results}
                  </div>
              </div>
          `;
      } else {
        search_result_container = `
            <div class="is-flex is-flex-direction-column gap-2 is-align-items-flex-start">
                ${modal_filters}
                ${search_divider}
                <div class="is-size-6">0 result(s)</div>
              </div>
              <div class="has-text-centered my-5 py-5">No result found!</div>
        `;
      }

      if ($(".search-modal-card-body").hasClass("is-justify-content-center")) {
        $(".search-modal-card-body").removeClass("is-justify-content-center");
      }

      $(".search-modal-card-body").html(search_result_container);
    } else {
      if (!$(".search-modal-card-body").hasClass("is-justify-content-center")) {
        $(".search-modal-card-body").addClass("is-justify-content-center");
      }

      $(".search-modal-card-body").html(`
        <div class="has-text-centered my-5 py-5">Type something to get started!</div>
      `);
    }
  }

  /**
   * Make the modal filter html
   *
   * @returns string
   */
  function make_modal_body_filters() {
    let str = filters
      .map((val) => {
        if (selected_filter == val.toLowerCase()) {
          return `<a href="javascript:;" class="search-filter search-filter-selected"><span>${val}</span></a>`;
        } else {
          return `<a href="javascript:;" class="search-filter"><span>${val}</span></a>`;
        }
      })
      .join("");

    return `
          <div class="is-flex gap-2 is-flex-wrap-wrap is-justify-content-flex-start is-align-items-center search-filters">
              <span class="is-size-6">Filters:</span>
              ${str}
          </div>`;
  }
}

function waitUntilSearchIndexAvailable() {
  // It is possible that the documenter.js script runs before the page
  // has finished loading and documenterSearchIndex gets defined.
  // So we need to wait until the search index actually loads before setting
  // up all the search-related stuff.
  if (typeof documenterSearchIndex !== "undefined") {
    runSearchMainCode();
  } else {
    console.warn("Search Index not available, waiting");
    setTimeout(waitUntilSearchIndexAvailable, 1000);
  }
}

// The actual entry point to the search code
waitUntilSearchIndexAvailable();

})
////////////////////////////////////////////////////////////////////////////////
require(['jquery'], function($) {

// Modal settings dialog
$(document).ready(function () {
  var settings = $("#documenter-settings");
  $("#documenter-settings-button").click(function () {
    settings.toggleClass("is-active");
  });
  // Close the dialog if X is clicked
  $("#documenter-settings button.delete").click(function () {
    settings.removeClass("is-active");
  });
  // Close dialog if ESC is pressed
  $(document).keyup(function (e) {
    if (e.keyCode == 27) settings.removeClass("is-active");
  });
});

})
////////////////////////////////////////////////////////////////////////////////
require(['jquery'], function($) {

$(document).ready(function () {
  let search_modal_header = `
    <header class="modal-card-head gap-2 is-align-items-center is-justify-content-space-between w-100 px-3">
      <div class="field mb-0 w-100">
        <p class="control has-icons-right">
          <input class="input documenter-search-input" type="text" placeholder="Search" />
          <span class="icon is-small is-right has-text-primary-dark">
            <i class="fas fa-magnifying-glass"></i>
          </span>
        </p>
      </div>
      <div class="icon is-size-4 is-clickable close-search-modal">
        <i class="fas fa-times"></i>
      </div>
    </header>
  `;

  let initial_search_body = `
    <div class="has-text-centered my-5 py-5">Type something to get started!</div>
  `;

  let search_modal_footer = `
    <footer class="modal-card-foot">
      <span>
        <kbd class="search-modal-key-hints">Ctrl</kbd> +
        <kbd class="search-modal-key-hints">/</kbd> to search
      </span>
      <span class="ml-3"> <kbd class="search-modal-key-hints">esc</kbd> to close </span>
    </footer>
  `;

  $(document.body).append(
    `
      <div class="modal" id="search-modal">
        <div class="modal-background"></div>
        <div class="modal-card search-min-width-50 search-min-height-100 is-justify-content-center">
          ${search_modal_header}
          <section class="modal-card-body is-flex is-flex-direction-column is-justify-content-center gap-4 search-modal-card-body">
            ${initial_search_body}
          </section>
          ${search_modal_footer}
        </div>
      </div>
    `
  );

  document.querySelector(".docs-search-query").addEventListener("click", () => {
    openModal();
  });

  document
    .querySelector(".close-search-modal")
    .addEventListener("click", () => {
      closeModal();
    });

  $(document).on("click", ".search-result-link", function () {
    closeModal();
  });

  document.addEventListener("keydown", (event) => {
    if ((event.ctrlKey || event.metaKey) && event.key === "/") {
      openModal();
    } else if (event.key === "Escape") {
      closeModal();
    }

    return false;
  });

  // Functions to open and close a modal
  function openModal() {
    let searchModal = document.querySelector("#search-modal");

    searchModal.classList.add("is-active");
    document.querySelector(".documenter-search-input").focus();
  }

  function closeModal() {
    let searchModal = document.querySelector("#search-modal");
    let initial_search_body = `
      <div class="has-text-centered my-5 py-5">Type something to get started!</div>
    `;

    searchModal.classList.remove("is-active");
    document.querySelector(".documenter-search-input").blur();

    if (!$(".search-modal-card-body").hasClass("is-justify-content-center")) {
      $(".search-modal-card-body").addClass("is-justify-content-center");
    }

    $(".documenter-search-input").val("");
    $(".search-modal-card-body").html(initial_search_body);
  }

  document
    .querySelector("#search-modal .modal-background")
    .addEventListener("click", () => {
      closeModal();
    });
});

})
////////////////////////////////////////////////////////////////////////////////
require(['jquery'], function($) {

// Manages the showing and hiding of the sidebar.
$(document).ready(function () {
  var sidebar = $("#documenter > .docs-sidebar");
  var sidebar_button = $("#documenter-sidebar-button");
  sidebar_button.click(function (ev) {
    ev.preventDefault();
    sidebar.toggleClass("visible");
    if (sidebar.hasClass("visible")) {
      // Makes sure that the current menu item is visible in the sidebar.
      $("#documenter .docs-menu a.is-active").focus();
    }
  });
  $("#documenter > .docs-main").bind("click", function (ev) {
    if ($(ev.target).is(sidebar_button)) {
      return;
    }
    if (sidebar.hasClass("visible")) {
      sidebar.removeClass("visible");
    }
  });
});

// Resizes the package name / sitename in the sidebar if it is too wide.
// Inspired by: https://github.com/davatron5000/FitText.js
$(document).ready(function () {
  e = $("#documenter .docs-autofit");
  function resize() {
    var L = parseInt(e.css("max-width"), 10);
    var L0 = e.width();
    if (L0 > L) {
      var h0 = parseInt(e.css("font-size"), 10);
      e.css("font-size", (L * h0) / L0);
      // TODO: make sure it survives resizes?
    }
  }
  // call once and then register events
  resize();
  $(window).resize(resize);
  $(window).on("orientationchange", resize);
});

// Scroll the navigation bar to the currently selected menu item
$(document).ready(function () {
  var sidebar = $("#documenter .docs-menu").get(0);
  var active = $("#documenter .docs-menu .is-active").get(0);
  if (typeof active !== "undefined") {
    sidebar.scrollTop = active.offsetTop - sidebar.offsetTop - 15;
  }
});

})
////////////////////////////////////////////////////////////////////////////////
require(['jquery'], function($) {

// Theme picker setup
$(document).ready(function () {
  // onchange callback
  $("#documenter-themepicker").change(function themepick_callback(ev) {
    var themename = $("#documenter-themepicker option:selected").attr("value");
    if (themename === "auto") {
      // set_theme(window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light');
      window.localStorage.removeItem("documenter-theme");
    } else {
      // set_theme(themename);
      window.localStorage.setItem("documenter-theme", themename);
    }
    // We re-use the global function from themeswap.js to actually do the swapping.
    set_theme_from_local_storage();
  });

  // Make sure that the themepicker displays the correct theme when the theme is retrieved
  // from localStorage
  if (typeof window.localStorage !== "undefined") {
    var theme = window.localStorage.getItem("documenter-theme");
    if (theme !== null) {
      $("#documenter-themepicker option").each(function (i, e) {
        e.selected = e.value === theme;
      });
    }
  }
});

})
////////////////////////////////////////////////////////////////////////////////
require(['jquery'], function($) {

// update the version selector with info from the siteinfo.js and ../versions.js files
$(document).ready(function () {
  // If the version selector is disabled with DOCUMENTER_VERSION_SELECTOR_DISABLED in the
  // siteinfo.js file, we just return immediately and not display the version selector.
  if (
    typeof DOCUMENTER_VERSION_SELECTOR_DISABLED === "boolean" &&
    DOCUMENTER_VERSION_SELECTOR_DISABLED
  ) {
    return;
  }

  var version_selector = $("#documenter .docs-version-selector");
  var version_selector_select = $("#documenter .docs-version-selector select");

  version_selector_select.change(function (x) {
    target_href = version_selector_select
      .children("option:selected")
      .get(0).value;
    window.location.href = target_href;
  });

  // add the current version to the selector based on siteinfo.js, but only if the selector is empty
  if (
    typeof DOCUMENTER_CURRENT_VERSION !== "undefined" &&
    $("#version-selector > option").length == 0
  ) {
    var option = $(
      "<option value='#' selected='selected'>" +
        DOCUMENTER_CURRENT_VERSION +
        "</option>"
    );
    version_selector_select.append(option);
  }

  if (typeof DOC_VERSIONS !== "undefined") {
    var existing_versions = version_selector_select.children("option");
    var existing_versions_texts = existing_versions.map(function (i, x) {
      return x.text;
    });
    DOC_VERSIONS.forEach(function (each) {
      var version_url = documenterBaseURL + "/../" + each + "/";
      var existing_id = $.inArray(each, existing_versions_texts);
      // if not already in the version selector, add it as a new option,
      // otherwise update the old option with the URL and enable it
      if (existing_id == -1) {
        var option = $(
          "<option value='" + version_url + "'>" + each + "</option>"
        );
        version_selector_select.append(option);
      } else {
        var option = existing_versions[existing_id];
        option.value = version_url;
        option.disabled = false;
      }
    });
  }

  // only show the version selector if the selector has been populated
  if (version_selector_select.children("option").length > 0) {
    version_selector.toggleClass("visible");
  }
});

})
