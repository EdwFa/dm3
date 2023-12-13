import React, { Component } from "react";
import { useState, useEffect, createRef } from "react";

import { Navigate, Link } from "react-router-dom";

import Graph from "react-graph-vis";
import { v4 as uuidv4 } from "uuid";
//import "./network.css";

import { AgGridReact } from "ag-grid-react";
import "ag-grid-enterprise";
import "ag-grid-community/styles/ag-grid.css";
import "ag-grid-community/styles/ag-theme-alpine.css";
import "../ag-theme-acmecorp.css";

import Plot from "react-plotly.js";

import { VOSviewerOnline } from "vosviewer-online";

import Select from "react-select";

import Slider from "react-input-slider";

import update from "immutability-helper";

import { variables, AG_GRID_LOCALE_RU } from "../Variables.js";

var topicFilterParams = {
  comparator: (TopicParam, cellValue) => {
    if (TopicParam === cellValue) {
      return 0;
    }
    if (cellValue < TopicParam) {
      return -1;
    }
    if (cellValue > TopicParam) {
      return -1;
    }
    return 0;
  },
};

const per_topics = [
  "Поиск в pubmed",
  "Тематический анализ",
  "Поиск в векторном представлении",
];

var ErrorMessage = 200;
var ErrorMessageText = "";
var Topic = "Все";

const obj_color = {
  disease: "#fdbbbb",
  drug: "#ECC58B",
  gene: "#E2DB8C",
  chemical: "#21c354",
  species: "#A6EFDC",
  mutation: "#B2DDEA",
  cell_type: "#C6DEF5",
  cell_line: "#A3B3D2",
  DNA: "#C9B9E8",
  RNA: "#D7DBE8",
};

function markup_text(text, annotations) {
  if (!annotations) {
    return text;
  }
  let markup_text = "";
  let last_position = 0;
  for (let annotation of annotations) {
    let start = annotation.span.begin;
    let end = annotation.span.end;
    if (!annotation.prop) {
      console.log("This");
    }
    let markup_str = `<span style=\"color: ${
      obj_color[annotation.obj]
    }\">${text.slice(start, end)}<sub>${
      annotation.prob ? annotation.prob.toFixed(2) : ""
    }</sub></span>`;
    markup_text = `${markup_text}${text.slice(
      last_position,
      start
    )}${markup_str}`;

    last_position = end;
  }
  return markup_text;
}

export class TematicReview extends Component {
  constructor(props) {
    super(props);

    this.gridRef = createRef();
    this.gridAnaliseRef = createRef();
    this.state = {
      updateOr: false,
      loading: false,
      useAll: true,
      token: variables.token,
      allow_page: variables.allow,

      // Search
      articles: [],
      DetailArticle: null,
      articlesInfo: [
        {
          field: "uid",
          filter: "agNumberColumnFilter",
          sortable: true,
          enableRowGroup: true,
          enableValue: true,
          resizable: true,
          headerName: "PMID",
        },
        {
          field: "titl",
          filter: "agTextColumnFilter",
          sortable: true,
          enableRowGroup: true,
          enableValue: true,
          minWidth: 300,
          width: 450,
          resizable: true,
          headerName: "Заголовок",
        },
        {
          field: "pdat",
          filter: "agTextColumnFilter",
          sortable: true,
          enableRowGroup: true,
          enableValue: true,
          resizable: true,
          headerName: "Дата выхода",
        },
        {
          field: "auth",
          filter: "agTextColumnFilter",
          sortable: true,
          enableRowGroup: true,
          enableValue: true,
          minWidth: 300,
          width: 450,
          resizable: true,
          headerName: "Авторы",
        },
        {
          field: "affl",
          filter: "agTextColumnFilter",
          sortable: true,
          enableRowGroup: true,
          enableValue: true,
          minWidth: 300,
          width: 450,
          resizable: true,
          headerName: "Аффилиации",
        },
        {
          field: "jour",
          filter: "agTextColumnFilter",
          sortable: true,
          enableRowGroup: true,
          enableValue: true,
          resizable: true,
          headerName: "Журнал",
        },
        {
          field: "pt",
          filter: "agTextColumnFilter",
          sortable: true,
          enableRowGroup: true,
          enableValue: true,
          resizable: true,
          headerName: "Тип статьи",
        },
        {
          field: "mesh",
          filter: "agTextColumnFilter",
          sortable: true,
          enableRowGroup: true,
          enableValue: true,
          minWidth: 300,
          width: 450,
          resizable: true,
          headerName: "MESH термины",
        },
      ],
      translation_stack: null,
      full_query: null,
      short_query: null,
      task: null,
      message: null,
      messageStatus: 200,
      count: 0,
      analiseRows: [],

      // Filters
      queryText: "",
      queryStartDate: "2022-01-01",
      queryEndDate: new Date().toISOString().split("T")[0],
      queryTypes: new Set(),
      queryOlds: new Set(),
      queryGenders: new Set(),

      // Analise
      messageAnalise: null,
      messageStatusAnalise: 200,

      // Analise filters
      rangeMin: 1,
      rangeMax: 3,
      min_TOPIC_SIZE: 10,
      top_N_WORDS: 10,

      top_n_topics: 40,
      n_clusters: 10,

      n_neighbors: 10,
      n_components: 2,
      min_dist: 0.0,
      metric: { label: "cosine" },
      list_of_metrics: [
        { label: "euclidean" },
        { label: "manhattan" },
        { label: "chebyshev" },
        { label: "minkowski" },
        { label: "canberra" },
        { label: "braycurtis" },
        { label: "mahalanobis" },
        { label: "wminkowski" },
        { label: "seuclidean" },
        { label: "cosine" },
        { label: "correlation" },
        { label: "haversine" },
        { label: "hamming" },
        { label: "jaccard" },
        { label: "dice" },
        { label: "russelrao" },
        { label: "kulsinski" },
        { label: "ll_dirichlet" },
        { label: "hellinger" },
        { label: "rogerstanimoto" },
        { label: "sokalmichener" },
        { label: "sokalsneath" },
        { label: "yule" },
      ],

      // Analise table
      analise_articles: [],
      analise_info: [
        {
          field: "titl",
          filter: "agTextColumnFilter",
          sortable: true,
          enableRowGroup: true,
          minWidth: 300,
          width: 450,
          resizable: true,
          headerName: "Заголовок",
        },
        {
          field: "pdat",
          filter: "agTextColumnFilter",
          sortable: true,
          enableRowGroup: true,
          resizable: true,
          headerName: "Дата выхода",
        },
        {
          field: "auth",
          filter: "agTextColumnFilter",
          sortable: true,
          enableRowGroup: true,
          minWidth: 300,
          width: 450,
          resizable: true,
          headerName: "Авторы",
        },
        {
          field: "affl",
          filter: "agTextColumnFilter",
          sortable: true,
          enableRowGroup: true,
          enableValue: true,
          minWidth: 300,
          width: 450,
          resizable: true,
          headerName: "Аффилиации",
        },
        {
          field: "jour",
          filter: "agTextColumnFilter",
          sortable: true,
          enableRowGroup: true,
          resizable: true,
          headerName: "Журнал",
        },
        {
          field: "pt",
          filter: "agTextColumnFilter",
          sortable: true,
          enableRowGroup: true,
          resizable: true,
          headerName: "Тип статьи",
        },
        {
          field: "mesh",
          filter: "agTextColumnFilter",
          sortable: true,
          enableRowGroup: true,
          minWidth: 300,
          width: 450,
          resizable: true,
          headerName: "MESH теримны",
        },
        {
          field: "topic",
          filter: "agNumberColumnFilter",
          sortable: true,
          enableRowGroup: true,
          filterParams: topicFilterParams,
          resizable: true,
          headerName: "Тема",
        },
        {
          field: "prop",
          filter: "agNumberColumnFilter",
          sortable: true,
          enableRowGroup: true,
          resizable: true,
          headerName: "Точность",
        },
      ],
      DetailArticle: null,

      // Analise clust graph
      clust_graph: null,
      heapmap: null,
      heirarchy: null,
      DTM: null,
      plotlyWidth: 800,

      // Filter topic
      current_topic: "Все",
      topicObject: { label: "Все" },
      topics: new Set(),

      //Summirise
      summarise: null,

      //Graphs
      messageGraph: null,
      messageStatusGraph: 200,

      current_graph: { label: "authors" },
      list_of_graphs: [
        { label: "authors" },
        { label: "affiliations" },
        { label: "journals" },
        { label: "countries" },
      ],

      infoAuthorsData: null,
      infoCountryData: null,
      infoJournalData: null,
      infoAffiliationsData: null,

      permissions: [],

      // Tabs
      useSearch: true,
      useAnalise: false,
      useGraphs: false,
    };
  }

  // Permissions
  getPermissions() {
    fetch(variables.API_URL + "/api/permissions", {
      headers: {
        "Content-Type": "application/json;charset=utf-8",
        Authorization: `Token ${variables.token}`,
      },
    })
      .then((response) => {
        console.log(response.status);
        ErrorMessage = response.status;
        if (response.ok) {
          return response.json();
        } else {
          throw Error(response.status);
        }
      })
      .then((data) => {
        console.log(data);
        this.setState({ permissions: data.permissions });
      })
      .catch((error) => {
        console.log(error);
      });
  }

  // Search

  getArticles = (url, interval = 1000) => {
    fetch(variables.API_URL + "/api/search/", {
      headers: {
        "Content-Type": "application/json;charset=utf-8",
        Authorization: `Token ${variables.token}`,
      },
    })
      .then((response) => {
        console.log(response.status);
        ErrorMessage = response.status;
        if (response.ok) {
          return response.json();
        } else {
          throw Error(response.status);
        }
      })
      .then((data) => {
        if (ErrorMessage === 202) {
          this.setState({
            loading: true,
            message: data.message,
            messageStatus: 202,
          });
          console.log(data.message);
          setTimeout(() => {
            return this.getArticles(url, interval);
          }, interval);
        } else {
          console.log("This work");
          console.log(
            data.task.query
              .split("AND ")[1]
              .split(":")[1]
              .replace("/", "-")
              .replace("/", "-")
          );
          this.setState({
            articles: data.data.search_ncbi,
            loading: false,
            full_query: data.task.full_query,
            translation_stack: data.task.translation_stack,
            short_query: data.task.query,
            count: data.task.count,
            message: "Запрос успешно обработан",
            messageStatus: 200,
            queryStartDate: data.task.query
              .split("AND ")[1]
              .split(":")[0]
              .replace("/", "-")
              .replace("/", "-"),
            queryEndDate: data.task.query
              .split("AND ")[1]
              .split(":")[1]
              .replace("/", "-")
              .replace("/", "-")
              .replace("[dp]", ""),
            queryText: data.task.query.split(" AND")[0],
          });
          this.getPermissions();
        }
      })
      .catch((error) => {
        if (ErrorMessage === 500) {
          this.setState({
            articles: [],
            DetailArticle: null,
            loading: false,
            message: "Ошибка сервера",
            messageStatus: 500,
          });
        } else if (ErrorMessage > 399) {
          this.setState({
            articles: [],
            DetailArticle: null,
            loading: false,
            message: "Что-то пошло не так A",
            messageStatus: 400,
          });
        }
      });
  };

  createTask() {
    // Отправляем запрос на сервер для получения статей
    if (this.state.queryText === "") {
      alert("Пожайлуста заполните поле запроса!");
      return;
    }
    this.setState({ loading: true });
    fetch(variables.API_URL + "/api/search/", {
      method: "POST",
      headers: {
        "Content-Type": "application/json;charset=utf-8",
        Authorization: `Token ${this.state.token}`,
      },
      body: JSON.stringify({
        search_field: this.state.queryText,
        dateStart: this.state.queryStartDate,
        dateStop: this.state.queryEndDate,
        Gender: [...this.state.queryGenders],
        Type: [...this.state.queryTypes],
        Old: [...this.state.queryOlds],
      }),
    })
      .then((response) => {
        console.log(response.status);
        if (response.ok) {
          return response.json();
        } else {
          ErrorMessage = response.status;
          throw Error(response.status);
        }
      })
      .then((data) => {
        this.setState({
          full_query: data.full_query,
          translation_stack: data.translation_stack,
          short_query: data.query,
          message: "Ваш запрос в очереди. Пожайлуста дождитесь результата",
          task: data,
          count: data.count,
          articles: [],
          loading: true,
          messageStatus: 201,
        });
        this.getArticles();
      })
      .catch((error) => {
        console.log(error);
        if (ErrorMessage === 500) {
          this.setState({
            task: null,
            loading: false,
            message: "Ошибка сервера",
            messageStatus: 500,
          });
        } else if (ErrorMessage === 403) {
          this.setState({
            task: null,
            loading: false,
            message: "Дождитесь окончания предыдушего запроса",
            messageStatus: 400,
          });
        } else {
          this.setState({
            task: null,
            loading: false,
            message: "Что-то пошло не так",
            messageStatus: 400,
          });
        }
      });
  }

  componentDidMount() {
    this.setState({
      message: "Запрашиваем данные с сервера",
      messageStatus: 202,
      loading: true,
      messageAnalise: "Запрашиваем данные с сервера",
      messageStatusAnalise: 202,
    });
    this.getArticles();
    this.getAnalise();
    this.getGraphInfo();
    console.log("start");
  }

  onSelectionAnalise = () => {
    const selectedRows = this.gridRef.current.api.getSelectedRows();
    this.setState({
      DetailArticle: selectedRows.length === 1 ? selectedRows[0] : null,
    });
  };

  changeQueryText = (e) => {
    this.setState({ queryText: e.target.value });
  };

  changeQueryStartDate = (e) => {
    this.setState({ queryStartDate: e.target.value });
  };

  changeQueryEndDate = (e) => {
    this.setState({ queryEndDate: e.target.value });
  };

  changeQueryGenders(gender) {
    if (this.state.queryGenders.has(gender)) {
      this.state.queryGenders.delete(gender);
    } else {
      this.state.queryGenders.add(gender);
    }
    this.setState({ updateOr: !this.state.updateOr });
  }

  changeQueryTypes(type) {
    if (this.state.queryTypes.has(type)) {
      this.state.queryTypes.delete(type);
    } else {
      this.state.queryTypes.add(type);
    }
    this.setState({ updateOr: !this.state.updateOr });
  }

  changeQueryOlds(old) {
    if (this.state.queryOlds.has(old)) {
      this.state.queryOlds.delete(old);
    } else {
      this.state.queryOlds.add(old);
    }
    this.setState({ updateOr: !this.state.updateOr });
  }

  changeQueryOldsMany(olds) {
    if (!this.state.useAll) {
      for (let old of olds) {
        this.state.queryOlds.delete(old);
      }
    } else {
      for (let old of olds) {
        this.state.queryOlds.add(old);
      }
    }
    this.setState({ useAll: !this.state.useAll });
  }

  changeModelRangeMin = (e) => {
    if (e.target.value > 5) {
      this.setState({ rangeMin: 5 });
    } else if (e.target.value < 1) {
      this.setState({ rangeMin: 1 });
    } else {
      this.setState({ rangeMin: e.target.value });
    }
  };

  changeModelRangeMax = (e) => {
    if (e.target.value > 5) {
      this.setState({ rangeMax: 5 });
    } else if (e.target.value < 1) {
      this.setState({ rangeMax: 1 });
    } else {
      this.setState({ rangeMax: e.target.value });
    }
  };

  changeMinDist = (e) => {
    if (e.target.value > 0.99) {
      this.setState({ min_dist: 0.99 });
    } else if (e.target.value < 0.0) {
      this.setState({ min_dist: 0.0 });
    } else {
      this.setState({ min_dist: e.target.value });
    }
  };

  RoundPersent(number) {
    return number.toFixed(2);
  }

  startAnalise() {
    let analise_data = [];
    this.gridRef.current.api.forEachNodeAfterFilter((rowNode) =>
      analise_data.push(rowNode.data.uid)
    );
    this.setState({ loading: true });
    fetch(variables.API_URL + "/api/analise/", {
      method: "POST",
      headers: {
        Accept: "application/json",
        "Content-Type": "application/json;charset=utf-8",
        Authorization: `Token ${variables.token}`,
      },
      body: JSON.stringify({
        articles: analise_data,
        filters: {
          rangeMin: this.state.rangeMin,
          rangeMax: this.state.rangeMax,
          min_TOPIC_SIZE: this.state.min_TOPIC_SIZE,
          top_N_WORDS: this.state.top_N_words,
          top_n_topics: this.state.top_n_topics,
          n_clusters: this.state.n_clusters,
          n_neighbors: this.state.n_neighbors,
          n_components: this.state.n_components,
          min_dist: this.state.min_dist,
          metric: this.state.metric.label,
        },
      }),
    })
      .then((response) => {
        console.log(response.status);
        if (response.ok) {
          return response.json();
        } else {
          ErrorMessage = response.status;
          throw Error(response.status);
        }
      })
      .then((data) => {
        this.setState({
          messageAnalise:
            "Ваш запрос в очереди. Пожайлуста дождитесь результата",
          messageStatusAnalise: 201,
          loading: true,
        });
        this.getAnalise();
      })
      .catch((error) => {
        if (ErrorMessage === 500) {
          this.setState({
            data: [],
            dataInfo: [],
            DetailArticle: null,
            loading: false,
            messageAnalise: "Ошибка сервера",
            messageStatusAnalise: 500,
          });
        } else if (ErrorMessage === 403) {
          this.setState({
            data: [],
            dataInfo: [],
            DetailArticle: null,
            loading: false,
            messageAnalise: "Дождитесь окончания предыдушего запроса",
            messageStatusAnalise: 403,
          });
        } else {
          this.setState({
            data: [],
            dataInfo: [],
            DetailArticle: null,
            loading: false,
            messageAnalise: "Вы пока не можете отправить запрос",
            messageStatusAnalise: 400,
          });
        }
      });
  }

  getAllArticles() {
    this.setState({loading: true});
    fetch(variables.API_URL + "/api/all_records/", {
      headers: {
        "Content-Type": "application/json;charset=utf-8",
        Authorization: `Token ${variables.token}`,
      },
    })
      .then((response) => {
        console.log(response.status);
        ErrorMessage = response.status;
        if (response.ok) {
          return response.json();
        } else {
          throw Error(response.status);
        }
      })
      .then((data) => {
          this.setState({
            articles: data.data.search_ncbi,
            loading: false,
          })
      })
      .catch((error) => {
        if (ErrorMessage === 500) {
          this.setState({
            articles: [],
            DetailArticle: null,
            loading: false,
            message: "Ошибка сервера",
            messageStatus: 500,
          });
        } else if (ErrorMessage > 399) {
          this.setState({
            articles: [],
            DetailArticle: null,
            loading: false,
            message: "Что-то пошло не так A",
            messageStatus: 400,
          });
        }
      });
  };

  // Analise

  onSelectionChanged = (gridApi) => {
    const selectedRows = this.gridAnaliseRef.current.api.getSelectedRows();
    this.setState({
      DetailArticle: selectedRows.length === 1 ? selectedRows[0] : null,
    });
  };

  getAnalise = (url, interval = 1000) => {
    fetch(variables.API_URL + "/api/analise/", {
      headers: {
        "Content-Type": "application/json;charset=utf-8",
        Authorization: `Token ${variables.token}`,
      },
    })
      .then((res) => {
        if (res.ok) {
          ErrorMessage = res.status;
          return res.json();
        } else {
          throw new Error(res);
        }
      })
      .then((data) => {
        console.log(data);
        if (ErrorMessage === 202) {
          this.setState({
            loading: true,
            messageAnalise: data.message,
            messageStatusAnalise: 202,
          });
          setTimeout(() => {
            return this.getAnalise(url, interval);
          }, interval);
        } else {
          if (data.data.clust_graph !== null) {
            // delete data.data.clust_graph.layout.width;
            data.data.clust_graph.layout.width = 800;
          }
          if (data.data.heapmap !== null)
            if (data.data.heapmap.length !== 0) {
              data.data.heapmap.layout.width = 800;
            }
          if (data.data.heirarchy !== null)
            if (data.data.heirarchy.length !== 0) {
              data.data.heirarchy.layout.width = 800;
            }
          if (data.data.DTM !== null)
            if (data.data.DTM.length !== 0) {
              data.data.DTM.layout.width = 800;
            }
          let topics = new Array();
          topics.push({ label: "Все" });
          for (let el of data.data.topics) {
            topics.push({ label: el });
          }
          console.log(topics);
          this.setState({
            analise_articles: data.data.tematic_analise,
            DetailArticle: data.data.tematic_analise[0],
            clust_graph: data.data.clust_graph,
            heapmap: data.data.heapmap,
            heirarchy: data.data.heirarchy,
            DTM: data.data.DTM,
            loading: false,
            messageAnalise: "Запрос успешно обработан",
            messageStatusAnalise: 200,
            topics: topics,
          });
          this.getPermissions();
        }
      })
      .catch((error) => {
        console.log(error);
        if (ErrorMessage === 500) {
          this.setState({
            analise_articles: [],
            clust_graph: null,
            heapmap: null,
            heirarchy: null,
            DetailArticle: null,
            loading: false,
            messageAnalise: "Ошибка сервера",
            messageStatusAnalise: 500,
          });
        } else if (ErrorMessage === 403) {
          this.setState({
            data: [],
            dataInfo: [],
            DetailArticle: null,
            loading: false,
            messageAnalise: "Дождитесь окончания предыдушего запроса",
            messageStatusAnalise: 403,
          });
        } else {
          this.setState({
            analise_articles: [],
            clust_graph: null,
            heapmap: null,
            heirarchy: null,
            DetailArticle: null,
            loading: false,
            messageAnalise: "Что-то пошло не так",
            messageStatusAnalise: 400,
          });
        }
      });
  };

  changePlotlyWidth(newWidth) {
    const heapmap = update(this.state.heapmap, {
      layout: { width: { $set: newWidth } },
    });

    const heirarchy = update(this.state.heirarchy, {
      layout: { width: { $set: newWidth } },
    });

    const clust_graph = update(this.state.clust_graph, {
      layout: { width: { $set: newWidth } },
    });

    const DTM = update(this.state.DTM, {
      layout: { width: { $set: newWidth } },
    });

    this.setState({
      plotlyWidth: newWidth,
      heapmap: heapmap,
      heirarchy: heirarchy,
      clust_graph: clust_graph,
      DTM: DTM,
    });
  }

  externalFilterChanged = (newValue) => {
    console.log(newValue);
    this.setState({
      current_topic: newValue.label.split("_")[0],
      summarise: null,
      topicObject: newValue,
    });
    if (newValue.label !== "Все") {
      Topic = Number(newValue.label.split("_")[0]);
    } else {
      Topic = newValue.label;
    }
    console.log(Topic);
    this.gridAnaliseRef.current.api.onFilterChanged();
  };

  isExternalFilterPresent = () => {
    // if ageType is not everyone, then we are filtering
    return Topic !== "Все";
  };

  doesExternalFilterPass = (node) => {
    if (node.data) {
      if (node.data.topic === Topic) {
        return true;
      }
    }
    return false;
  };

  autoGroupColumnDef = () => {
    return {
      minWidth: 200,
    };
  };

  getGraphData = () => {
    var current_topic = this.state.current_topic;
    if (current_topic === "Все") {
      return this.state.clust_graph.data;
    }

    var data = [];
    for (let topic of this.state.clust_graph.data) {
      var topic_id = topic.name.split("_")[0];
      if (current_topic === topic_id) {
        data.push(topic);
      }
    }
    return data;
  };

  getSummarise = (task_id, interval = 1000) => {
    fetch(variables.API_URL + `/api/summarise?task_id=${task_id}`, {
      headers: {
        "Content-Type": "application/json;charset=utf-8",
        Authorization: `Token ${variables.token}`,
      },
    })
      .then((res) => {
        if (res.status == 202) {
          this.setState({
            loading: true,
            messageStatusAnalise: 202,
            messageAnalise: "Отправлено на суммаризацию...",
          });
          setTimeout(() => {
            return this.getSummarise(task_id, interval);
          }, interval);
        }
        if (res.status == 200) {
          return res.json();
        } else {
          ErrorMessage = res.status;
          throw Error(res.statusText);
        }
      })
      .then((data) => {
        this.setState({
          summarise: data.data,
          loading: false,
          messageAnalise: "Суммаризация прошла успешно",
          messageStatusAnalise: 200,
        });
      })
      .catch((err) => {
        console.log(err);
        if (ErrorMessage === 202) {
          this.setState({
            loading: true,
            messageStatusAnalise: 202,
            messageAnalise: "Отправлено на суммаризацию...",
          });
        } else {
          this.setState({
            loading: false,
            summarise: null,
            messageAnalise: "Произошла ошибка при суммаризации",
            messageStatusAnalise: 500,
          });
        }
      });
  };

  createSummariseQuery() {
    let data = [];
    this.gridAnaliseRef.current.api.forEachNodeAfterFilter((rowNode) =>
      data.push(rowNode.data.uid)
    );
    this.setState({ loading: true });
    fetch(variables.API_URL + "/api/summarise", {
      method: "POST",
      headers: {
        Accept: "application/json",
        "Content-Type": "application/json;charset=utf-8",
        Authorization: `Token ${variables.token}`,
      },
      body: JSON.stringify({
        articles: data,
      }),
    })
      .then((response) => {
        console.log(response.status);
        if (response.ok) {
          return response.json();
        } else {
          ErrorMessage = response.status;
          throw Error(response.status);
        }
      })
      .then((result) => {
        var task_id = result.data;
        this.setState({
          messageAnalise:
            "Отправлено на суммаризацию пожайлуста дождитесь ответа",
          messageStatusAnalise: 201,
          loading: true,
        });
        this.getSummarise(task_id);
      })
      .catch((error) => {
        this.setState({
          messageAnalise: "Ошибка при суммаризации",
          messageStatusAnalise: 500,
          loading: false,
        });
      });
  }

  // Graphs authors, countries, jornals

  getGraphInfo = (url, interval = 1000) => {
    fetch(variables.API_URL + "/api/graphs/", {
      headers: {
        "Content-Type": "application/json;charset=utf-8",
        Authorization: `Token ${variables.token}`,
      },
    })
      .then((response) => {
        console.log(response.status);
        ErrorMessage = response.status;
        if (response.ok) {
          return response.json();
        } else {
          throw Error(response.status);
        }
      })
      .then((data) => {
        if (ErrorMessage === 202) {
          this.setState({
            loading: true,
            messageAnalise: data.message,
            messageStatusAnalise: 202,
          });
          console.log(data.message);
          setTimeout(() => {
            return this.getGraphInfo(url, interval);
          }, interval);
        } else {
          console.log("This graph work");
          this.setState({
            loading: false,
            infoAuthorsData: data.data.info_graph,
            infoAffiliationsData: data.data.info_graph_affiliations,
            infoJournalData: data.data.info_graph_journals,
            infoCountryData: data.data.info_graph_countries,
            messageAnalise:
              "Граф успешно отрисован, перейдите во вкладку графы для просмотра",
            messageStatusAnalise: 200,
          });
        }
      })
      .catch((error) => {
        if (ErrorMessage === 500) {
          this.setState({
            articles: [],
            DetailArticle: null,
            loading: false,
            messageAnalise: "Ошибка сервера",
            messageStatusAnalise: 500,
          });
        } else {
          this.setState({
            articles: [],
            DetailArticle: null,
            loading: false,
            messageAnalise: "Что-то пошло не так",
            messageStatusAnalise: 400,
          });
        }
      });
  };

  createGraph() {
    let analise_data = [];
    this.gridAnaliseRef.current.api.forEachNodeAfterFilter((rowNode) =>
      analise_data.push(rowNode.data.uid)
    );
    this.setState({ loading: true });
    fetch(variables.API_URL + "/api/graphs/", {
      method: "POST",
      headers: {
        Accept: "application/json",
        "Content-Type": "application/json;charset=utf-8",
        Authorization: `Token ${variables.token}`,
      },
      body: JSON.stringify({
        articles: analise_data,
      }),
    })
      .then((response) => {
        console.log(response.status);
        if (response.ok) {
          return response.json();
        } else {
          ErrorMessage = response.status;
          throw Error(response.status);
        }
      })
      .then((data) => {
        this.setState({
          messageAnalise:
            "Ваш запрос в очереди. Пожайлуста дождитесь результата",
          messageStatusAnalise: 201,
          loading: true,
        });
        this.getGraphInfo();
      })
      .catch((error) => {
        if (ErrorMessage === 500) {
          this.setState({
            data: [],
            dataInfo: [],
            DetailArticle: null,
            loading: false,
            messageAnalise: "Ошибка сервера",
            messageStatusAnalise: 500,
          });
        } else if (ErrorMessage === 403) {
          this.setState({
            data: [],
            dataInfo: [],
            DetailArticle: null,
            loading: false,
            messageAnalise: "Дождитесь окончания предыдушего запроса",
            messageStatusAnalise: 403,
          });
        } else {
          this.setState({
            data: [],
            dataInfo: [],
            DetailArticle: null,
            loading: false,
            messageAnalise: "Что=то пошло не так",
            messageStatusAnalise: 400,
          });
        }
      });
  }

  // разметка

  getMarkUp = (task_id, interval = 1000) => {
    fetch(variables.API_URL + `/api/markup?task_id=${task_id}`, {
      headers: {
        "Content-Type": "application/json;charset=utf-8",
        Authorization: `Token ${variables.token}`,
      },
    })
      .then((res) => {
        if (res.status == 202) {
          this.setState({ loading: true });
          setTimeout(() => {
            return this.getMarkUp(task_id, interval);
          }, interval);
        } else if (res.status == 200) {
          return res.json();
        } else {
          throw Error(res.statusText);
        }
      })
      .then((data) => {
        try {
          this.setState({
            DetailArticle: data.data,
            message: "Разметка прошла успешно",
            messageStatus: 200,
            loading: false,
          });
        } catch {
          console.log("access");
        }
      })
      .catch((err) => {
        console.log(err);
        this.setState({
          message: "Произошла ошибка при разметке",
          loading: false,
          messageStatus: 500,
        });
      });
  };

  markUpArticle(DetailArticle) {
    this.setState({ loading: true });
    fetch(variables.API_URL + "/api/markup", {
      method: "POST",
      headers: {
        Accept: "application/json",
        "Content-Type": "application/json;charset=utf-8",
        Authorization: `Token ${variables.token}`,
      },
      body: JSON.stringify({
        article: DetailArticle,
      }),
    })
      .then((res) => {
        console.log(res.status);
        if (res.ok) {
          return res.json();
        } else {
          throw Error(res.statusText);
        }
      })
      .then((result) => {
        var task_id = result.data;
        this.setState({
          message: "Отправлено на суммаризацию пожайлуста дождитесь ответа",
          messageStatus: 201,
          loading: true,
        });
        this.getMarkUp(task_id);
      })
      .catch((err) => {
        console.log(err);
        this.setState({
          message: "ошибка при разметке",
          messageStatus: 500,
          loading: false,
        });
      });
  }

  // Скачивание embeddings
  downloadEmbeddings() {
    fetch(variables.API_URL + "/api/download_vectors", {
      method: "GET",
      headers: {
        Authorization: `Token ${this.state.token}`,
      },
      responseType: "blob",
    })
      .then((res) => {
        if (res.status == 200) {
          console.log(res);
          res.blob().then((blob) => {
            let url = window.URL.createObjectURL(blob);
            let a = document.createElement("a");
            a.href = url;
            a.download = "vectors_OR.tsv";
            a.click();
          });
        } else {
          throw Error(res.statusText);
        }
      })
      .catch((error) => {
        alert("Ошибка");
        console.log(error);
      });
  }

  downloadMetadata() {
    fetch(variables.API_URL + "/api/download_metadata", {
      method: "GET",
      headers: {
        Authorization: `Token ${this.state.token}`,
      },
      responseType: "blob",
    })
      .then((res) => {
        if (res.status == 200) {
          console.log(res);
          res.blob().then((blob) => {
            let url = window.URL.createObjectURL(blob);
            let a = document.createElement("a");
            a.href = url;
            a.download = "metadata_OR.tsv";
            a.click();
          });
        } else {
          throw Error(res.statusText);
        }
      })
      .catch((error) => {
        alert("Ошибка");
        console.log(error);
      });
  }

  downloadAll() {
    this.downloadMetadata();
    this.downloadEmbeddings();
  }

  render() {
    const {
      token,
      count,
      articles,
      DetailArticle,
      articlesInfo,
      translation_stack,
      full_query,
      short_query,
      message,
      messageStatus,
      messageAnalise,
      messageStatusAnalise,
      loading,

      queryText,
      queryEndDate,
      queryStartDate,

      rangeMin,
      rangeMax,
      top_n_topics,
      n_components,
      n_neighbors,
      n_clusters,
      min_TOPIC_SIZE,
      top_N_WORDS,
      min_dist,
      metric,
      list_of_metrics,

      analise_articles,
      analise_info,
      clust_graph,
      heapmap,
      heirarchy,
      DTM,
      current_topic,
      topics,
      topicObject,
      summarise,
      allow_page,

      current_graph,
      list_of_graphs,
      infoAuthorsData,
      infoCountryData,
      infoJournalData,
      infoAffiliationsData,

      messageGraph,
      messageStatusGraph,
      plotlyWidth,

      permissions,

      useSearch,
      useAnalise,
      useGraphs,
    } = this.state;

    if (!token) {
      return <Navigate push to="/login" />;
    } else if (allow_page === 1) {
      return <Navigate push to="/ddi_review" />
    } else if (allow_page === 2) {
      return <Navigate push to="/chat" />
    } else {
      return (
        <>
          <header className="bg-white">
            <nav class="px-2 py-2.5">
              <div class="flex flex-wrap justify-between items-center">
                <div class="flex justify-start items-center">
                  <a href="" class="flex mr-4">
                    <img
                      src="https://flowbite.s3.amazonaws.com/logo.svg"
                      class="mr-3 h-8"
                      alt="FlowBite Logo"
                    />
                    <span class="self-center text-2xl font-semibold whitespace-nowrap">
                      EBM Sechenov DataMed.AI
                    </span>
                  </a>
                  {allow_page === 3 ? (
                    <ul class="flex font-medium flex-row space-x-8 ml-10">
                      <Link to="/tematic_review">
                        <li>
                          <a
                            href="#"
                            class="block py-2 pl-3 pr-4 text-gray-900 bg-blue-700 rounded md:bg-transparent md:text-blue-700 md:p-0"
                            aria-current="page"
                          >
                            Тематический анализ
                          </a>
                        </li>
                      </Link>
                      <Link to="/ddi_review">
                        <li>
                          <a
                            href="#"
                            class="block py-2 pl-3 pr-4 text-gray-900 rounded hover:bg-gray-100 md:hover:bg-transparent md:hover:text-blue-700 md:p-0"
                          >
                            Факты для EBM
                          </a>
                        </li>
                      </Link>
                      {variables.admin?
                        <Link to="/admin">
                          <li>
                            <a href="#" class="block py-2 pl-3 pr-4 text-gray-900 rounded hover:bg-gray-100 md:hover:bg-transparent md:hover:text-blue-700 md:p-0">Админ панель</a>
                          </li>
                        </Link>
                      :null}
                    </ul>
                  ) : null}
                </div>
                <div class="flex items-center lg:order-2">
                  <div class="flex-shrink-0 dropdown">
                    <a
                      href="#"
                      class="d-block link-body-emphasis text-decoration-none dropdown-toggle"
                      data-bs-toggle="dropdown"
                      aria-expanded="false"
                    >
                      <img
                        src="https://github.com/mdo.png"
                        alt="mdo"
                        width="32"
                        height="32"
                        class="rounded-circle"
                      />
                    </a>
                    <ul class="dropdown-menu text-small shadow">
                      {permissions?.map((per) => (
                        <li>
                          <a class="dropdown-item" href="#">
                            {per.topic}{" "}
                            {per.all_records
                              ? `${per.all_records}`
                              : "безлимитно"}
                          </a>
                        </li>
                      ))}
                    </ul>
                  </div>
                </div>
              </div>
            </nav>
            <nav class="px-2 py-2 border border-gray-200">
              <div class="w-full">
                <div class="flex justify-between items-center">
                  <button
                    id="toggleSidebar"
                    aria-expanded="true"
                    aria-controls="sidebar"
                    class="hidden order-first p-2 mr-3 text-gray-600 rounded cursor-pointer lg:inline hover:text-gray-900 hover:bg-gray-100 dark:text-gray-400 dark:hover:text-white dark:hover:bg-gray-700"
                    data-bs-toggle="collapse"
                    data-bs-target="#sidebar"
                    aria-label="Toggle navigation"
                  >
                    <svg
                      class="w-6 h-6"
                      fill="currentColor"
                      viewBox="0 0 20 20"
                      xmlns="http://www.w3.org/2000/svg"
                    >
                      <path
                        fill-rule="evenodd"
                        d="M3 5a1 1 0 011-1h12a1 1 0 110 2H4a1 1 0 01-1-1zM3 10a1 1 0 011-1h6a1 1 0 110 2H4a1 1 0 01-1-1zM3 15a1 1 0 011-1h12a1 1 0 110 2H4a1 1 0 01-1-1z"
                        clip-rule="evenodd"
                      ></path>
                    </svg>
                  </button>
                  <div className="grow">
                    <div class="relative">
                      <div class="absolute inset-y-0 left-0 flex items-center pl-3 pointer-events-none">
                        <svg
                          class="w-4 h-4 text-gray-500 dark:text-gray-400"
                          aria-hidden="true"
                          xmlns="http://www.w3.org/2000/svg"
                          fill="none"
                          viewBox="0 0 20 20"
                        >
                          <path
                            stroke="currentColor"
                            stroke-linecap="round"
                            stroke-linejoin="round"
                            stroke-width="2"
                            d="m19 19-4-4m0-7A7 7 0 1 1 1 8a7 7 0 0 1 14 0Z"
                          />
                        </svg>
                      </div>
                      <input
                        class="w-full py-3 bg-gray-50 border border-gray-300 text-gray-900 sm:text-sm rounded-lg focus:ring-primary-500 focus:border-primary-500 block w-full pl-10 p-2.5"
                        id="search"
                        type="text"
                        name="search_field"
                        placeholder={
                          short_query
                            ? short_query.split(" AND")[0]
                            : "covid-19"
                        }
                        value={queryText}
                        onChange={this.changeQueryText}
                        aria-label="Search"
                      />
                      <button
                        type="submit"
                        value="Найти"
                        disabled={loading}
                        onClick={() => this.createTask()}
                        class="text-white absolute right-2.5 bottom-2.5 bg-blue-700 hover:bg-blue-800 focus:ring-4 focus:outline-none focus:ring-blue-300 font-medium rounded-lg text-sm px-4 py-2"
                      >
                        {loading?
                                        <svg
                                          class="-ml-1 ml-3 h-5 w-5 text-white"
                                          xmlns="http://www.w3.org/2000/svg"
                                          fill="none"
                                          viewBox="0 0 24 24"
                                        >
                                          <circle
                                            class="opacity-25"
                                            cx="12"
                                            cy="12"
                                            r="10"
                                            stroke="currentColor"
                                            stroke-width="4"
                                          ></circle>
                                          <path
                                            class="opacity-75"
                                            fill="currentColor"
                                            d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"
                                          ></path>
                                        </svg>
                                :
                                    'Найти'
                                }
                      </button>
                    </div>
                  </div>
                  <div
                    class="ml-5 items-center justify-between hidden w-full md:flex md:w-auto md:order-1"
                    id="navbar-sticky"
                  >
                    <ul class="nav nav-pills" id="myTab" role="tablist">
                      <li class="nav-item mr-2" role="presentation">
                        <button
                          class="nav-link inline-block px-4 py-2 rounded-lg hover:text-gray-900 hover:bg-gray-100 active"
                          id="home-tab"
                          data-bs-toggle="tab"
                          data-bs-target="#home"
                          type="button"
                          role="tab"
                          aria-controls="home"
                          aria-selected={useSearch}
                          onClick={() =>
                            this.setState({
                              useGraphs: false,
                              useSearch: true,
                              useAnalise: false,
                            })
                          }
                        >
                          Результаты поиска
                        </button>
                      </li>
                      <li class="nav-item mr-2" role="presentation">
                        <button
                          class="nav-link inline-block px-4 py-2 rounded-lg hover:text-gray-900 hover:bg-gray-100"
                          id="profile-tab"
                          data-bs-toggle="tab"
                          data-bs-target="#profile"
                          type="button"
                          role="tab"
                          aria-controls="profile"
                          aria-selected={useAnalise}
                          onClick={() =>
                            this.setState({
                              useGraphs: false,
                              useSearch: false,
                              useAnalise: true,
                            })
                          }
                        >
                          Тематическое описание коллекции
                        </button>
                      </li>
                      <li class="nav-item mr-2" role="presentation">
                        <button
                          class="nav-link inline-block px-4 py-2 rounded-lg hover:text-gray-900 hover:bg-gray-100"
                          id="contact-tab"
                          data-bs-toggle="tab"
                          data-bs-target="#contact"
                          type="button"
                          role="tab"
                          aria-controls="contact"
                          aria-selected={useGraphs}
                          onClick={() =>
                            this.setState({
                              useGraphs: true,
                              useSearch: false,
                              useAnalise: false,
                            })
                          }
                        >
                          Графы аффиляций
                        </button>
                      </li>
                    </ul>
                    <button
                      id="toggleSidebar"
                      aria-expanded="true"
                      aria-controls="sidebar2"
                      class="order-last hidden p-2 text-gray-600 rounded cursor-pointer lg:inline hover:text-gray-900 hover:bg-gray-100 dark:text-gray-400 dark:hover:text-white dark:hover:bg-gray-700"
                      data-bs-toggle="collapse"
                      data-bs-target="#sidebar2"
                      aria-label="Toggle navigation"
                    >
                      <svg
                        class="w-6 h-6 rotate-180"
                        fill="currentColor"
                        viewBox="0 0 20 20"
                        xmlns="http://www.w3.org/2000/svg"
                      >
                        <path
                          fill-rule="evenodd"
                          d="M3 5a1 1 0 011-1h12a1 1 0 110 2H4a1 1 0 01-1-1zM3 10a1 1 0 011-1h6a1 1 0 110 2H4a1 1 0 01-1-1zM3 15a1 1 0 011-1h12a1 1 0 110 2H4a1 1 0 01-1-1z"
                          clip-rule="evenodd"
                        ></path>
                      </svg>
                    </button>
                  </div>
                </div>
              </div>
            </nav>
          </header>
          <main>
            <div>
              <div className="container-fluid h-screen">
                <div className="row align-items-stretch b-height">
                  <aside
                    id="sidebar"
                    className="h-screen col-md-2 my-3 bg-white collapse show width border rounded-3 overflow-auto g-0"
                  >
                    <div
                      className="accordion accordion-flush"
                      id="accordionFlushExample"
                    >
                      <div className="accordion-item">
                        <h2 className="accordion-header" id="flush-headingOne">
                          <button
                            className="accordion-button collapsed"
                            type="button"
                            data-bs-toggle="collapse"
                            data-bs-target="#flush-collapseOne"
                            aria-expanded="false"
                            aria-controls="flush-collapseOne"
                          >
                            Дата публикации
                          </button>
                        </h2>
                        <div
                          id="flush-collapseOne"
                          className="collapse show multi-collapse"
                          aria-labelledby="flush-headingOne"
                          data-bs-target="#accordionFlushExample"
                        >
                          <div className="accordion-body">
                            <div className="mb-3">
                              <label for="localdate">От : </label>
                              <input
                                type="date"
                                id="d1"
                                name="dateStart"
                                value={queryStartDate}
                                onChange={this.changeQueryStartDate}
                              />
                            </div>
                            <div>
                              <label for="localdate">До : </label>
                              <input
                                type="date"
                                id="d2"
                                name="dateStop"
                                value={queryEndDate}
                                onChange={this.changeQueryEndDate}
                              />
                            </div>
                          </div>
                        </div>
                      </div>
                      <div className="accordion-item">
                        <h2 class="accordion-header" id="flush-headingThree">
                          <button
                            class="accordion-button collapsed"
                            data-target="#flush-collapseThree"
                            type="button"
                            data-bs-toggle="collapse"
                            data-bs-target="#flush-collapseFour"
                            aria-expanded="false"
                            aria-controls="flush-collapseThree"
                          >
                            Тип статьи
                          </button>
                        </h2>
                        <div
                          id="flush-collapseFour"
                          className="collapse show multi-collapse"
                          aria-labelledby="flush-headingFour"
                          data-bs-target="#accordionFlushExample"
                        >
                          <div className="accordion-body">
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxClinicalTrial"
                                name="CheckBoxClinicalTrial"
                                checked={this.state.queryTypes.has(
                                  "clinical trial"
                                )}
                                onChange={() =>
                                  this.changeQueryTypes("clinical trial")
                                }
                              />
                              <label
                                className="form-check-label"
                                for="CheckboxClinicalTrial"
                              >
                                Clinical Trial
                              </label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxMetaAnalysys"
                                name="CheckboxMetaAnalysys"
                                checked={this.state.queryTypes.has(
                                  "meta-analysis"
                                )}
                                onChange={() =>
                                  this.changeQueryTypes("meta-analysis")
                                }
                              />
                              <label
                                className="form-check-label"
                                for="CheckboxMetaAnalysys"
                              >
                                Meta Analysys
                              </label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxRandomizedControlledTrial"
                                name="CheckboxRandomizedControlledTrial"
                                checked={this.state.queryTypes.has(
                                  "randomized controlled trial"
                                )}
                                onChange={() =>
                                  this.changeQueryTypes(
                                    "randomized controlled trial"
                                  )
                                }
                              />
                              <label
                                className="form-check-label"
                                for="CheckboxRandomizedControlledTrial"
                              >
                                Randomized Controlled Trial
                              </label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxReview"
                                name="CheckboxReview"
                                checked={this.state.queryTypes.has("review")}
                                onChange={() => this.changeQueryTypes("review")}
                              />
                              <label
                                className="form-check-label"
                                for="CheckboxReview"
                              >
                                Review
                              </label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxSystematicReview"
                                name="CheckboxSystematicReview"
                                checked={this.state.queryTypes.has(
                                  "systematic review"
                                )}
                                onChange={() =>
                                  this.changeQueryTypes("systematic review")
                                }
                              />
                              <label
                                className="form-check-label"
                                for="CheckboxSystematicReview"
                              >
                                Systematic Review
                              </label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxJournalArticle"
                                name="CheckboxJournalArticle"
                                checked={this.state.queryTypes.has(
                                  "journal article"
                                )}
                                onChange={() =>
                                  this.changeQueryTypes("journal article")
                                }
                              />
                              <label
                                className="form-check-label"
                                for="CheckboxJournalArticle"
                              >
                                Journal Article
                              </label>
                            </div>
                          </div>
                        </div>
                      </div>
                      <div className="accordion-item">
                        <h2
                          className="accordion-header"
                          id="flush-headingThree"
                        >
                          <button
                            className="accordion-button collapsed"
                            data-target="#flush-collapseThree"
                            type="button"
                            data-bs-toggle="collapse"
                            data-bs-target="#flush-collapseThree"
                            aria-expanded="false"
                            aria-controls="flush-collapseThree"
                          >
                            Возраст пациента
                          </button>
                        </h2>
                        <div
                          id="flush-collapseThree"
                          className="collapse show multi-collapse"
                          aria-labelledby="flush-headingThree"
                          data-bs-target="#accordionFlushExample"
                        >
                          <div className="accordion-body">
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxChild18"
                                checked={
                                  this.state.queryOlds.has("infant[mh]") &
                                  this.state.queryOlds.has("child[mh]") &
                                  this.state.queryOlds.has("adolescent[mh]")
                                }
                                onChange={() =>
                                  this.changeQueryOldsMany([
                                    "infant[mh]",
                                    "child[mh]",
                                    "adolescent[mh]",
                                  ])
                                }
                              />
                              <label
                                className="form-check-label"
                                for="CheckboxChild18"
                              >
                                Child: birth-18 years
                              </label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxNewborn"
                                checked={this.state.queryOlds.has(
                                  "infant, newborn[mh]"
                                )}
                                onChange={() =>
                                  this.changeQueryOlds("infant, newborn[mh]")
                                }
                              />
                              <label
                                className="form-check-label"
                                for="CheckboxNewborn"
                              >
                                Newborn: birth-1 months
                              </label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxInfant023"
                                checked={this.state.queryOlds.has("infant[mh]")}
                                onChange={() =>
                                  this.changeQueryOlds("infant[mh]")
                                }
                              />
                              <label
                                className="form-check-label"
                                for="CheckboxInfant023"
                              >
                                Infant: birth-23 months
                              </label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxInfant123"
                                checked={this.state.queryOlds.has(
                                  "infant[mh:noexp]"
                                )}
                                onChange={() =>
                                  this.changeQueryOlds("infant[mh:noexp]")
                                }
                              />
                              <label
                                className="form-check-label"
                                for="CheckboxInfant123"
                              >
                                Infant: 1-23 months
                              </label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxPreschool"
                                checked={this.state.queryOlds.has(
                                  "child, preschool[mh]"
                                )}
                                onChange={() =>
                                  this.changeQueryOlds("child, preschool[mh]")
                                }
                              />
                              <label
                                className="form-check-label"
                                for="CheckboxPreschool"
                              >
                                Preschool Child: 2-5 years
                              </label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxChild612"
                                checked={this.state.queryOlds.has(
                                  "child[mh:noexp]"
                                )}
                                onChange={() =>
                                  this.changeQueryOlds("child[mh:noexp]")
                                }
                              />
                              <label
                                className="form-check-label"
                                for="CheckboxChild612"
                              >
                                Child: 6-12 years
                              </label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxAdolescent"
                                checked={this.state.queryOlds.has(
                                  "adolescent[mh]"
                                )}
                                onChange={() =>
                                  this.changeQueryOlds("adolescent[mh]")
                                }
                              />
                              <label
                                className="form-check-label"
                                for="CheckboxAdolescent"
                              >
                                Adolescent: 13-18 years
                              </label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxAdult19"
                                checked={this.state.queryOlds.has("adult[mh]")}
                                onChange={() =>
                                  this.changeQueryOlds("adult[mh]")
                                }
                              />
                              <label
                                className="form-check-label"
                                for="CheckboxAdult19"
                              >
                                Adult: 19+ years
                              </label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxYoungAdult"
                                checked={this.state.queryOlds.has(
                                  "young adult[mh]"
                                )}
                                onChange={() =>
                                  this.changeQueryOlds("young adult[mh]")
                                }
                              />
                              <label
                                className="form-check-label"
                                for="CheckboxYoungAdult"
                              >
                                Young Adult: 19-24 years
                              </label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxAdult1944"
                                checked={this.state.queryOlds.has(
                                  "adult[mh:noexp]"
                                )}
                                onChange={() =>
                                  this.changeQueryOlds("adult[mh:noexp]")
                                }
                              />
                              <label
                                className="form-check-label"
                                for="CheckboxAdult1944"
                              >
                                Adult: 19-44 years
                              </label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxMiddle45"
                                checked={
                                  this.state.queryOlds.has("middle aged[mh]") &
                                  this.state.queryOlds.has("aged[mh]")
                                }
                                onChange={() =>
                                  this.changeQueryOldsMany([
                                    "middle aged[mh]",
                                    "aged[mh]",
                                  ])
                                }
                              />
                              <label
                                className="form-check-label"
                                for="CheckboxMiddle45"
                              >
                                Middle Aged: + Aged 45+ years
                              </label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxMiddle4565"
                                checked={this.state.queryOlds.has(
                                  "middle aged[mh]"
                                )}
                                onChange={() =>
                                  this.changeQueryOlds("middle aged[mh]")
                                }
                              />
                              <label
                                className="form-check-label"
                                for="CheckboxMiddle4565"
                              >
                                Middle Aged: 45-65 years
                              </label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxAged"
                                checked={this.state.queryOlds.has("aged[mh]")}
                                onChange={() =>
                                  this.changeQueryOlds("aged[mh]")
                                }
                              />
                              <label
                                className="form-check-label"
                                for="CheckboxAged"
                              >
                                Aged: 65+ years
                              </label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="Checkbox80"
                                checked={this.state.queryOlds.has(
                                  "aged, 80 and over[mh]"
                                )}
                                onChange={() =>
                                  this.changeQueryOlds("aged, 80 and over[mh]")
                                }
                              />
                              <label
                                className="form-check-label"
                                for="Checkbox80"
                              >
                                80 and over: 80+ years
                              </label>
                            </div>
                          </div>
                        </div>
                      </div>
                      <div className="accordion-item">
                        <h2 className="accordion-header" id="flush-headingTwo">
                          <button
                            className="accordion-button collapsed"
                            type="button"
                            data-bs-toggle="collapse"
                            data-bs-target="#flush-collapseTwo"
                            aria-expanded="false"
                            aria-controls="flush-collapseTwo"
                          >
                            Пол
                          </button>
                        </h2>
                        <div
                          id="flush-collapseTwo"
                          className="collapse show multi-collapse"
                          aria-labelledby="flush-headingTwo"
                          data-bs-target="#accordionFlushExample"
                        >
                          <div className="accordion-body">
                            <div className="form-check form-check-inline">
                              <input
                                type="checkbox"
                                id="CheckboxMan"
                                name="CheckboxMan"
                                className="form-check-input"
                                checked={this.state.queryGenders.has("man")}
                                onChange={() => this.changeQueryGenders("man")}
                              />
                              <label
                                className="form-check-label"
                                for="CheckboxMan"
                              >
                                Мужчина
                              </label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                name="CheckboxWomen"
                                id="CheckboxWomen"
                                checked={this.state.queryGenders.has("woman")}
                                onChange={() =>
                                  this.changeQueryGenders("woman")
                                }
                              />
                              <label
                                className="form-check-label"
                                for="CheckboxWomen"
                              >
                                Женщина
                              </label>
                            </div>
                          </div>
                        </div>
                      </div>
                      <div className="accordion-item">
                        <h2 className="accordion-header" id="flush-headingOne">
                          <button
                            className="accordion-button collapsed"
                            type="button"
                            data-bs-toggle="collapse"
                            data-bs-target="#flush-collapseTwelve"
                            aria-expanded="false"
                            aria-controls="flush-collapseOne"
                          >
                            Настройки Модели
                          </button>
                        </h2>
                        <div
                          id="flush-collapseTwelve"
                          className="collapse show multi-collapse"
                          aria-labelledby="flush-headingOne"
                          data-bs-target="#accordionFlushExample"
                        >
                          <div className="accordion-body">
                            <div className="form-outline">
                              <label class="form-label" for="typeNumber1">
                                Range Min:{" "}
                              </label>
                              <input
                                min="1"
                                max="5"
                                type="number"
                                id="typeNumber1"
                                className="form-control"
                                value={rangeMin}
                                onChange={this.changeModelRangeMin}
                              />
                            </div>
                            <div classNams="form-outline">
                              <label class="form-label" for="typeNumber2">
                                Range Max:{" "}
                              </label>
                              <input
                                min="1"
                                max="5"
                                type="number"
                                id="typeNumber2"
                                className="form-control"
                                value={rangeMax}
                                onChange={this.changeModelRangeMax}
                              />
                            </div>
                            <div classNams="form-outline">
                              <label class="form-label" for="typeNumber3">
                                min_TOPIC_SIZE:
                              </label>
                              <input
                                type="number"
                                id="typeNumber3"
                                className="form-control"
                                value={min_TOPIC_SIZE}
                                onChange={(x) =>
                                  this.setState({
                                    min_TOPIC_SIZE: x.target.value,
                                  })
                                }
                              />
                            </div>
                            <div classNams="form-outline">
                              <label class="form-label" for="typeNumber4">
                                top_N_WORDS:{" "}
                              </label>
                              <input
                                type="number"
                                id="typeNumber4"
                                className="form-control"
                                value={top_N_WORDS}
                                onChange={(x) =>
                                  this.setState({ top_N_WORDS: x.target.value })
                                }
                              />
                            </div>
                          </div>
                        </div>
                      </div>
                      <div className="accordion-item">
                        <h2 className="accordion-header" id="flush-headingOne">
                          <button
                            className="accordion-button collapsed"
                            type="button"
                            data-bs-toggle="collapse"
                            data-bs-target="#flush-collapseThirteen"
                            aria-expanded="false"
                            aria-controls="flush-collapseOne"
                          >
                            Матрица связаности
                          </button>
                        </h2>
                        <div
                          id="flush-collapseThirteen"
                          className="collapse show multi-collapse"
                          aria-labelledby="flush-headingOne"
                          data-bs-target="#accordionFlushExample"
                        >
                          <div className="accordion-body">
                            <div classNams="form-outline">
                              <label class="form-label" for="typeNumber5">
                                top_n_topics:{" "}
                              </label>
                              <input
                                type="number"
                                id="typeNumber5"
                                className="form-control"
                                value={top_n_topics}
                                onChange={(x) =>
                                  this.setState({
                                    top_n_topics: x.target.value,
                                  })
                                }
                              />
                            </div>
                            <div classNams="form-outline">
                              <label class="form-label" for="typeNumber6">
                                n_clusters:{" "}
                              </label>
                              <input
                                type="number"
                                id="typeNumber6"
                                className="form-control"
                                value={n_clusters}
                                onChange={(x) =>
                                  this.setState({ n_clusters: x.target.value })
                                }
                              />
                            </div>
                          </div>
                        </div>
                      </div>
                      <div className="accordion-item">
                        <h2 className="accordion-header" id="flush-headingOne">
                          <button
                            className="accordion-button collapsed"
                            type="button"
                            data-bs-toggle="collapse"
                            data-bs-target="#flush-collapseFourteen"
                            aria-expanded="false"
                            aria-controls="flush-collapseOne"
                          >
                            UMAP группировка
                          </button>
                        </h2>
                        <div
                          id="flush-collapseFourteen"
                          className="collapse show multi-collapse"
                          aria-labelledby="flush-headingOne"
                          data-bs-target="#accordionFlushExample"
                        >
                          <div className="accordion-body">
                            <div classNams="form-outline">
                              <label class="form-label" for="typeNumber7">
                                n_neighbors:{" "}
                              </label>
                              <input
                                type="number"
                                id="typeNumber7"
                                className="form-control"
                                value={n_neighbors}
                                onChange={(x) =>
                                  this.setState({ n_neighbors: x.target.value })
                                }
                              />
                            </div>
                            <div classNams="form-outline">
                              <label class="form-label" for="typeNumber8">
                                n_components:{" "}
                              </label>
                              <input
                                type="number"
                                id="typeNumber8"
                                className="form-control"
                                value={n_components}
                                onChange={(x) =>
                                  this.setState({
                                    n_components: x.target.value,
                                  })
                                }
                              />
                            </div>
                            <div classNams="form-outline">
                              <label class="form-label" for="typeNumber9">
                                min_dist:{" "}
                              </label>
                              <input
                                type="number"
                                id="typeNumber9"
                                className="form-control"
                                value={min_dist}
                                onChange={this.changeMinDist}
                              />
                            </div>
                            <div classNams="form-outline">
                              <label class="form-label" for="typeNumber9">
                                metric:{" "}
                              </label>
                              <Select
                                className="basic-single"
                                classNamePrefix="select"
                                value={metric}
                                isSearchable
                                placeholder="Выберите класс"
                                name="topic"
                                options={list_of_metrics}
                                getOptionLabel={(option) => option.label}
                                getOptionValue={(option) => option.label}
                                onChange={(x) => this.setState({ metric: x })}
                              />
                            </div>
                          </div>
                        </div>
                      </div>
                    </div>
                  </aside>
                  <section class="col p-3 m-3 border rounded-3 bg-white h-screen overflow-auto">
                    <div class="bd-example">
                      <div class="tab-content" id="myTabContent">
                        <div
                          class="tab-pane fade active show"
                          id="home"
                          role="tabpanel"
                          aria-labelledby="home-tab"
                        >
                          <div class="container-fluid g-0">
                            <div
                              class="accordion accordion-flush"
                              id="accordion"
                            >
                              <div class="accordion-item">
                                <div className="flex flex-row">
                                  <div
                                    class="accordion-header"
                                    id=""
                                    className="p-2 mr-2 rounded-lg bg-gray-200 border border-gray-500 grow font-medium text-sm "
                                  >
                                    {message ? (
                                      messageStatus > 299 ? (
                                        <p
                                          class="animate-pulse"
                                          style={{ color: "red" }}
                                        >
                                          {message}.{" "}
                                        </p>
                                      ) : messageStatus === 200 ? (
                                        <p
                                          class="animate-pulse"
                                          style={{ color: "green" }}
                                        >
                                          {message}.
                                        </p>
                                      ) : (
                                        <p
                                          class="animate-pulse"
                                          style={{ color: "black" }}
                                        >
                                          {message}.
                                        </p>
                                      )
                                    ) : null}
                                  </div>
                                  <button
                                    type="button"
                                    class="inline-flex items-center px-4 py-2 font-semibold leading-6 text-sm shadow rounded-md text-white bg-blue-700 hover:bg-blue-800 transition ease-in-out duration-150"
                                    disabled={loading}
                                    onClick={() => this.startAnalise()}
                                  >
                                    {loading?
                                        <svg
                                          class="-ml-1 ml-3 h-5 w-5 text-white"
                                          xmlns="http://www.w3.org/2000/svg"
                                          fill="none"
                                          viewBox="0 0 24 24"
                                        >
                                          <circle
                                            class="opacity-25"
                                            cx="12"
                                            cy="12"
                                            r="10"
                                            stroke="currentColor"
                                            stroke-width="4"
                                          ></circle>
                                          <path
                                            class="opacity-75"
                                            fill="currentColor"
                                            d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"
                                          ></path>
                                        </svg>
                                    :
                                        'Обработать'
                                    }
                                  </button>
                                </div>

                                <button
                                  class="accordion-button collapsed mt-2"
                                  type="button"
                                  data-bs-toggle="collapse"
                                  data-bs-target="#flush-collapseEleven"
                                  aria-expanded="false"
                                  aria-controls="flush-collapseSeven"
                                >
                                  По запросу найдено {count} источников.
                                </button>
                                <div
                                  id="flush-collapseEleven"
                                  class="collapse multi-collapse"
                                  aria-labelledby="flush-headingEleven"
                                  data-bs-target="#accordionFlushExample"
                                >
                                  <div class="accordion-body">
                                    <p class="pb-2 mb-3 border-bottom">
                                      {" "}
                                      Запрос {short_query} .
                                    </p>
                                    <p class="pb-2 mb-3 border-bottom">
                                      {" "}
                                      Запрос автоматически расширен до
                                      следующего вида - {full_query}.
                                    </p>
                                    <p class="pb-2 mb-3 border-bottom">
                                      {" "}
                                      Служебная информация для анализа :{" "}
                                      {translation_stack}.
                                    </p>
                                  </div>
                                </div>
                              </div>
                            </div>

                            <div
                              className="ag-theme-alpine ag-theme-acmecorp"
                              style={{ height: 700 }}
                            >
                              <AgGridReact
                                ref={this.gridRef}
                                rowData={articles}
                                columnDefs={articlesInfo}
                                pagination={true}
                                onSelectionChanged={this.onSelectionAnalise}
                                autoGroupColumnDef={this.autoGroupColumnDef}
                                onChange={this.externalFilterChanged}
                                rowSelection={"single"}
                                localeText={AG_GRID_LOCALE_RU}
                                sideBar={{
                                  toolPanels: [
                                    {
                                      id: "columns",
                                      labelDefault: "Columns",
                                      labelKey: "columns",
                                      iconKey: "columns",
                                      toolPanel: "agColumnsToolPanel",
                                      minWidth: 225,
                                      width: 225,
                                      maxWidth: 225,
                                    },
                                    {
                                      id: "filters",
                                      labelDefault: "Filters",
                                      labelKey: "filters",
                                      iconKey: "filter",
                                      toolPanel: "agFiltersToolPanel",
                                      minWidth: 180,
                                      maxWidth: 400,
                                      width: 250,
                                    },
                                  ],
                                  position: "left",
                                }}
                              ></AgGridReact>
                              <br />
                              <button
                                    type="button"
                                    class="inline-flex items-center px-4 py-2 font-semibold leading-6 text-sm shadow rounded-md text-white bg-blue-700 hover:bg-blue-800 transition ease-in-out duration-150"
                                    disabled={loading}
                                    onClick={() => this.getAllArticles()}
                              >
                                {loading?
                                        <svg
                                          class="-ml-1 ml-3 h-5 w-5 text-white"
                                          xmlns="http://www.w3.org/2000/svg"
                                          fill="none"
                                          viewBox="0 0 24 24"
                                        >
                                          <circle
                                            class="opacity-25"
                                            cx="12"
                                            cy="12"
                                            r="10"
                                            stroke="currentColor"
                                            stroke-width="4"
                                          ></circle>
                                          <path
                                            class="opacity-75"
                                            fill="currentColor"
                                            d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"
                                          ></path>
                                        </svg>
                                :
                                    'Загрузить полностью'
                                }
                              </button>

                            </div>
                          </div>
                        </div>
                        <div
                          class="tab-pane fade"
                          id="profile"
                          role="tabpanel"
                          aria-labelledby="profile-tab"
                        >
                          <div>
                            <h2 class="accordion-header" id="">
                              {messageAnalise ? (
                                messageStatusAnalise > 299 ? (
                                  <p
                                    class="pb-2 mb-3 border-bottom"
                                    style={{ color: "red" }}
                                  >
                                    {messageAnalise}.{" "}
                                  </p>
                                ) : messageStatusAnalise === 200 ? (
                                  <p
                                    class="pb-2 mb-3 border-bottom"
                                    style={{ color: "green" }}
                                  >
                                    {messageAnalise}.
                                  </p>
                                ) : (
                                  <p
                                    class="pb-2 mb-3 border-bottom"
                                    style={{ color: "black" }}
                                  >
                                    {messageAnalise}.
                                  </p>
                                )
                              ) : null}
                            </h2>
                            <Select
                              className="basic-single"
                              classNamePrefix="select"
                              value={topicObject}
                              isSearchable
                              placeholder="Выберите класс"
                              name="topic"
                              options={topics}
                              getOptionLabel={(option) => option.label}
                              getOptionValue={(option) => option.label}
                              onChange={this.externalFilterChanged}
                            />
                            <br />
                            <button
                                    type="button"
                                    className="text-white right-2.5 my-4 bottom-2.5 bg-blue-700 hover:bg-blue-800 focus:ring-4 focus:outline-none focus:ring-blue-300 font-medium rounded-lg text-sm px-4 py-2"
                                    disabled={loading}
                                    onClick={() => this.createGraph()}
                              >
                                {loading?
                                        <svg
                                          class="-ml-1 ml-3 h-5 w-5 text-white"
                                          xmlns="http://www.w3.org/2000/svg"
                                          fill="none"
                                          viewBox="0 0 24 24"
                                        >
                                          <circle
                                            class="opacity-25"
                                            cx="12"
                                            cy="12"
                                            r="10"
                                            stroke="currentColor"
                                            stroke-width="4"
                                          ></circle>
                                          <path
                                            class="opacity-75"
                                            fill="currentColor"
                                            d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"
                                          ></path>
                                        </svg>
                                :
                                    'Отрисовать граф'
                                }
                              </button>
                          </div>
                          <div class="accordion accordion-flush" id="accordion">
                            <div className="accordion-item">
                              <h2
                                className="accordion-header"
                                id="flush-headingTwo"
                              >
                                <button
                                  className="accordion-button collapsed"
                                  type="button"
                                  data-bs-toggle="collapse"
                                  data-bs-target="#flush-collapseSix"
                                  aria-expanded="false"
                                  aria-controls="flush-collapseTwo"
                                >
                                  Таблица
                                </button>
                              </h2>
                              <div
                                id="flush-collapseSix"
                                className="collapse show multi-collapse"
                                aria-labelledby="flush-headingTwo"
                                data-bs-target="#accordionFlushExample"
                              >
                                <div className="accordion-body">
                                  <div
                                    className="ag-theme-alpine ag-theme-acmecorp"
                                    style={{ height: 700 }}
                                  >
                                    <AgGridReact
                                      ref={this.gridAnaliseRef}
                                      rowData={analise_articles}
                                      localeText={AG_GRID_LOCALE_RU}
                                      columnDefs={analise_info}
                                      pagination={true}
                                      rowSelection={"single"}
                                      onSelectionChanged={
                                        this.onSelectionChanged
                                      }
                                      animateRows={true}
                                      isExternalFilterPresent={
                                        this.isExternalFilterPresent
                                      }
                                      doesExternalFilterPass={
                                        this.doesExternalFilterPass
                                      }
                                      sideBar={{
                                        toolPanels: [
                                          {
                                            id: "columns",
                                            labelDefault: "Columns",
                                            labelKey: "columns",
                                            iconKey: "columns",
                                            toolPanel: "agColumnsToolPanel",
                                            minWidth: 225,
                                            width: 225,
                                            maxWidth: 225,
                                          },
                                          {
                                            id: "filters",
                                            labelDefault: "Filters",
                                            labelKey: "filters",
                                            iconKey: "filter",
                                            toolPanel: "agFiltersToolPanel",
                                            minWidth: 180,
                                            maxWidth: 400,
                                            width: 250,
                                          },
                                        ],
                                        position: "left",
                                      }}
                                    ></AgGridReact>
                                  </div>
                                </div>
                              </div>
                            </div>
                          </div>
                          <label class="form-label">
                            Ширина = {plotlyWidth}
                          </label>
                          <Slider
                            axis="x"
                            x={plotlyWidth}
                            xmax={2000}
                            xmin={200}
                            xstep={100}
                            onChange={({ x }) => this.changePlotlyWidth(x)}
                          />
                          <div class="accordion accordion-flush" id="accordion">
                            <div className="accordion-item">
                              <h2
                                className="accordion-header"
                                id="flush-headingTwo"
                              >
                                <button
                                  className="accordion-button collapsed"
                                  type="button"
                                  data-bs-toggle="collapse"
                                  data-bs-target="#flush-collapseSeven"
                                  aria-expanded="false"
                                  aria-controls="flush-collapseTwo"
                                >
                                  Тематические коллекции
                                </button>
                              </h2>
                              <div
                                id="flush-collapseSeven"
                                className="collapse show multi-collapse"
                                aria-labelledby="flush-headingTwo"
                                data-bs-target="#accordionFlushExample"
                              >
                                <div className="accordion-body">
                                  <div>
                                    {clust_graph ? (
                                      <Plot
                                        data={this.getGraphData()}
                                        layout={clust_graph.layout}
                                        style={{
                                          width: "100%",
                                          height: "100%",
                                        }}
                                        useResizeHandler={true}
                                      />
                                    ) : null}
                                  </div>
                                </div>
                              </div>
                            </div>

                            <div className="accordion-item">
                              <h2
                                className="accordion-header"
                                id="flush-headingTwo"
                              >
                                <button
                                  className="accordion-button collapsed"
                                  type="button"
                                  data-bs-toggle="collapse"
                                  data-bs-target="#flush-collapseEight"
                                  aria-expanded="false"
                                  aria-controls="flush-collapseTwo"
                                >
                                  Матрица близости тем
                                </button>
                              </h2>
                              <div
                                id="flush-collapseEight"
                                className="collapse show multi-collapse"
                                aria-labelledby="flush-headingTwo"
                                data-bs-target="#accordionFlushExample"
                              >
                                <div className="accordion-body">
                                  <div>
                                    {heapmap ? (
                                      <Plot
                                        data={heapmap.data}
                                        layout={heapmap.layout}
                                      />
                                    ) : null}
                                  </div>
                                </div>
                              </div>
                            </div>
                            <div className="accordion-item">
                              <h2
                                className="accordion-header"
                                id="flush-headingTwo"
                              >
                                <button
                                  className="accordion-button collapsed"
                                  type="button"
                                  data-bs-toggle="collapse"
                                  data-bs-target="#flush-collapseNine"
                                  aria-expanded="false"
                                  aria-controls="flush-collapseTwo"
                                >
                                  Иерархия тем
                                </button>
                              </h2>
                              <div
                                id="flush-collapseNine"
                                className="collapse show multi-collapse"
                                aria-labelledby="flush-headingTwo"
                                data-bs-target="#accordionFlushExample"
                              >
                                <div className="accordion-body">
                                  <div>
                                    {heirarchy ? (
                                      <Plot
                                        data={heirarchy.data}
                                        layout={heirarchy.layout}
                                      />
                                    ) : null}
                                  </div>
                                </div>
                              </div>
                            </div>
                            <div className="accordion-item">
                              <h2
                                className="accordion-header"
                                id="flush-headingTwo"
                              >
                                <button
                                  className="accordion-button collapsed"
                                  type="button"
                                  data-bs-toggle="collapse"
                                  data-bs-target="#flush-collapseTen"
                                  aria-expanded="false"
                                  aria-controls="flush-collapseTwo"
                                >
                                  Временная динамика по темам
                                </button>
                              </h2>
                              <div
                                id="flush-collapseTen"
                                className="collapse show multi-collapse"
                                aria-labelledby="flush-headingTwo"
                                data-bs-target="#accordionFlushExample"
                              >
                                <div className="accordion-body">
                                  <div>
                                    {DTM ? (
                                      <Plot
                                        data={DTM.data}
                                        layout={DTM.layout}
                                      />
                                    ) : null}
                                  </div>
                                </div>
                              </div>
                            </div>
                            <div className="accordion-item">
                              <h2
                                className="accordion-header"
                                id="flush-headingTwo"
                              >
                                <button
                                  className="accordion-button collapsed"
                                  type="button"
                                  data-bs-toggle="collapse"
                                  data-bs-target="#flush-collapseEleven"
                                  aria-expanded="false"
                                  aria-controls="flush-collapseTwo"
                                >
                                  Суммаризация
                                </button>

                              </h2>
                              <div
                                id="flush-collapseEleven"
                                className="collapse show multi-collapse"
                                aria-labelledby="flush-headingTwo"
                                data-bs-target="#accordionFlushExample"
                              >
                                <div className="accordion-body">
                                  <div>
                                    {summarise ? (
                                      <>
                                        <p>Summarise</p>
                                        <p>{summarise}</p>
                                      </>
                                    ) : (
                                      <button
                                        type="button"
                                        className="text-white right-2.5 my-4 bottom-2.5 bg-blue-700 hover:bg-blue-800 focus:ring-4 focus:outline-none focus:ring-blue-300 font-medium rounded-lg text-sm px-4 py-2"
                                        disabled={loading}
                                        onClick={() => this.createSummariseQuery()}
                                        >
                                        {loading?
                                                <svg
                                                  class="-ml-1 ml-3 h-5 w-5 text-white"
                                                  xmlns="http://www.w3.org/2000/svg"
                                                  fill="none"
                                                  viewBox="0 0 24 24"
                                                >
                                                  <circle
                                                    class="opacity-25"
                                                    cx="12"
                                                    cy="12"
                                                    r="10"
                                                    stroke="currentColor"
                                                    stroke-width="4"
                                                  ></circle>
                                                  <path
                                                    class="opacity-75"
                                                    fill="currentColor"
                                                    d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"
                                                  ></path>
                                                </svg>
                                        :
                                            'Суммаризовать'
                                        }
                                        </button>
                                    )}
                                  </div>
                                </div>
                              </div>
                            </div>
                            <div className="accordion-item">
                              <h2
                                className="accordion-header"
                                id="flush-headingTwo"
                              >
                                <button
                                  className="accordion-button collapsed"
                                  type="button"
                                  data-bs-toggle="collapse"
                                  data-bs-target="#flush-collapseZero"
                                  aria-expanded="false"
                                  aria-controls="flush-collapseTwo"
                                >
                                  Проектор
                                </button>
                              </h2>
                              <div
                                id="flush-collapseZero"
                                className="collapse show multi-collapse"
                                aria-labelledby="flush-headingTwo"
                                data-bs-target="#accordionFlushExample"
                              >
                                <div className="accordion-body">
                                  <div>
                                    <a
                                      href="http://projector.tensorflow.org/"
                                      class="card-title link-primary text-decoration-none h6"
                                      target="_blank"
                                    >
                                      {" "}
                                      Перейти на проектор{" "}
                                    </a>
                                    <button
                                        type="button"
                                        className="text-white right-2.5 my-4 bottom-2.5 bg-blue-700 hover:bg-blue-800 focus:ring-4 focus:outline-none focus:ring-blue-300 font-medium rounded-lg text-sm px-4 py-2"
                                        disabled={loading}
                                        onClick={() => this.downloadAll()}
                                    >
                                    {loading?
                                            <svg
                                              class="-ml-1 ml-3 h-5 w-5 text-white"
                                              xmlns="http://www.w3.org/2000/svg"
                                              fill="none"
                                              viewBox="0 0 24 24"
                                            >
                                              <circle
                                                class="opacity-25"
                                                cx="12"
                                                cy="12"
                                                r="10"
                                                stroke="currentColor"
                                                stroke-width="4"
                                              ></circle>
                                              <path
                                                class="opacity-75"
                                                fill="currentColor"
                                                d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"
                                              ></path>
                                            </svg>
                                    :
                                        'Cкачать'
                                    }
                                    </button>
                                  </div>
                                </div>
                              </div>
                            </div>
                          </div>
                        </div>
                        <div
                          class="tab-pane fade"
                          id="contact"
                          role="tabpanel"
                          aria-labelledby="contact-tab"
                        >
                          <div class="container-fluid g-0">
                            <h2 class="accordion-header" id="">
                              {messageGraph ? (
                                messageStatusGraph > 299 ? (
                                  <p
                                    class="pb-2 mb-3 border-bottom"
                                    style={{ color: "red" }}
                                  >
                                    {messageGraph}.{" "}
                                  </p>
                                ) : messageStatusGraph === 200 ? (
                                  <p
                                    class="pb-2 mb-3 border-bottom"
                                    style={{ color: "green" }}
                                  >
                                    {messageGraph}.
                                  </p>
                                ) : (
                                  <p
                                    class="pb-2 mb-3 border-bottom"
                                    style={{ color: "black" }}
                                  >
                                    {messageGraph}.
                                  </p>
                                )
                              ) : null}
                            </h2>
                            <Select
                              className="basic-single"
                              classNamePrefix="select"
                              value={current_graph}
                              isSearchable
                              placeholder="Выберите класс"
                              name="topic"
                              options={list_of_graphs}
                              getOptionLabel={(option) => option.label}
                              getOptionValue={(option) => option.label}
                              onChange={(x) =>
                                this.setState({ current_graph: x })
                              }
                            />
                            <div
                              id="mynetwork"
                              style={{ width: "100%", height: "600px" }}
                            >
                              {current_graph.label === "authors" ? (
                                infoAuthorsData ? (
                                  <VOSviewerOnline data={infoAuthorsData} />
                                ) : null
                              ) : null}
                              {current_graph.label === "affiliations" ? (
                                infoAffiliationsData ? (
                                  <VOSviewerOnline
                                    data={infoAffiliationsData}
                                  />
                                ) : null
                              ) : null}
                              {current_graph.label === "countries" ? (
                                infoCountryData ? (
                                  <VOSviewerOnline data={infoCountryData} />
                                ) : null
                              ) : null}
                              {current_graph.label === "journals" ? (
                                infoCountryData ? (
                                  <VOSviewerOnline data={infoJournalData} />
                                ) : null
                              ) : null}
                            </div>
                          </div>
                        </div>
                      </div>
                    </div>
                  </section>

                  <aside
                    id="sidebar2"
                    class="col-md-4 h-screen collapse show width col p-3 my-3 border rounded-3 overflow-auto bg-white"
                  >
                    <h3 class="pb-2 mb-3 border-bottom">Подробное описание</h3>
                    <nav class="small" id="toc">
                      {DetailArticle ? (
                        <div class="card mb-3">
                          <div class="card-body">
                            <a
                              href={DetailArticle.url}
                              class="card-title link-primary text-decoration-none h5"
                              target="_blank"
                            >
                              {" "}
                              {DetailArticle.titl}{" "}
                            </a>
                            <p class="card-text">
                              ----------------------------------{" "}
                            </p>
                            <p class="card-text">
                              Авторы : {DetailArticle.auth}{" "}
                            </p>
                            <p class="card-text">
                              ----------------------------------{" "}
                            </p>
                            <p class="card-text">Аннотация : </p>
                            <p
                              class="card-text"
                              dangerouslySetInnerHTML={{
                                __html: markup_text(
                                  DetailArticle.tiab,
                                  DetailArticle.annotations
                                ),
                              }}
                            />
                            <p class="card-text">
                              ----------------------------------{" "}
                            </p>
                            <p class="card-text">
                              <small class="text-success">
                                Дата публикации : {DetailArticle.pdat}{" "}
                              </small>
                            </p>
                            <p class="card-text">
                              <small class="text-success">
                                Издание : {DetailArticle.jour}
                              </small>
                            </p>
                            <p class="card-text">
                              <small class="text-success">
                                Вид публикации : {DetailArticle.pt}
                              </small>
                            </p>
                            <p class="card-text">
                              <small class="text-success">
                                Страна : {DetailArticle.pl}{" "}
                              </small>
                            </p>
                            <p class="card-text">
                              <small class="text-success">
                                {DetailArticle.mesh}{" "}
                              </small>
                            </p>
                            {loading ? (
                              <p>Loading...</p>
                            ) : (
                              <input
                                className="text-white right-2.5 my-4 bottom-2.5 bg-blue-700 hover:bg-blue-800 focus:ring-4 focus:outline-none focus:ring-blue-300 font-medium rounded-lg text-sm px-4 py-2"
                                type="submit"
                                value="Разметить"
                                disabled={loading}
                                onClick={() =>
                                  this.markUpArticle(DetailArticle)
                                }
                              />
                            )}
                          </div>
                        </div>
                      ) : null}
                    </nav>
                  </aside>
                </div>
              </div>
            </div>
          </main>
        </>
      );
    }
  }
}