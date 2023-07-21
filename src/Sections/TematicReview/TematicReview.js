import React, { Component } from 'react';
import { useState, useEffect, createRef } from 'react';

import { Navigate, Link } from 'react-router-dom';

import Graph from "react-graph-vis";
import { v4 as uuidv4 } from 'uuid'
//import "./network.css";

import { AgGridReact } from 'ag-grid-react';
import 'ag-grid-enterprise';
import 'ag-grid-community/styles/ag-grid.css';
import 'ag-grid-community/styles/ag-theme-alpine.css';
import '../ag-theme-acmecorp.css';

import Plot from 'react-plotly.js';

import { VOSviewerOnline } from 'vosviewer-online'

import Select from 'react-select';

import Slider from 'react-input-slider';

import { variables } from '../Variables.js';



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

var ErrorMessage = 0

const obj_color = {
  'disease': "#fdbbbb",
  'drug': "#ECC58B",
  'gene': "#E2DB8C",
  'chemical': "#21c354",
  'species': "#A6EFDC",
  'mutation': "#B2DDEA",
  'cell_type': "#C6DEF5",
  'cell_line': "#A3B3D2",
  'DNA': "#C9B9E8",
  'RNA': "#D7DBE8",
}

function markup_text(text, annotations) {
  if (!annotations) {
    return text
  }
  let markup_text = ''
  let last_position = 0
  for (let annotation of annotations) {
    let start = annotation.span.begin
    let end = annotation.span.end
    if (!annotation.prop) { console.log('This') }
    let markup_str = `<span style=\"color: ${obj_color[annotation.obj]}\">${text.slice(start, end)}<sub>${annotation.prob ? annotation.prob.toFixed(2) : ""}</sub></span>`
    markup_text = `${markup_text}${text.slice(last_position, start)}${markup_str}`

    last_position = end
  }
  return markup_text
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
        { field: 'titl', filter: 'agTextColumnFilter', enableValue: true, minWidth: 300, width: 450, resizable: true},
        { field: 'pdat', filter: 'agTextColumnFilter', enableRowGroup: true, enableValue: true, resizable: true},
        { field: 'auth', filter: 'agTextColumnFilter', enableValue: true, minWidth: 300, width: 450, resizable: true},
        { field: 'jour', filter: 'agTextColumnFilter', enableRowGroup: true, enableValue: true, resizable: true},
        { field: 'pt', filter: 'agTextColumnFilter', enableRowGroup: true, enableValue: true, resizable: true},
        { field: 'mesh', filter: 'agTextColumnFilter', enableValue: true, minWidth: 300, width: 450, resizable: true},
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
      queryText: '',
      queryStartDate: '2022-01-01',
      queryEndDate: new Date().toISOString().split('T')[0],
      queryTypes: new Set(),
      queryOlds: new Set(),
      queryGenders: new Set(),

      // Analise
      messageAnalise: null,
      messageStatusAnalise: 200,

      // Analise table
      analise_articles: [],
      analise_info: [
        { field: 'titl', filter: 'agTextColumnFilter', minWidth: 300, width: 450, resizable: true},
        { field: 'pdat', filter: 'agTextColumnFilter', resizable: true},
        { field: 'auth', filter: 'agTextColumnFilter', minWidth: 300, width: 450, resizable: true},
        { field: 'jour', filter: 'agTextColumnFilter', resizable: true},
        { field: 'pt', filter: 'agTextColumnFilter', resizable: true},
        { field: 'mesh', filter: 'agTextColumnFilter', minWidth: 300, width: 450, resizable: true},
        { field: 'topic', filter: 'agNumberColumnFilter', sortable: true, filterParams: topicFilterParams, resizable: true},
        { field: 'prop', filter: 'agNumberColumnFilter', resizable: true},
      ],
      DetailArticle: null,

      // Analise clust graph
      clust_graph: null,
      heapmap: null,
      heirarchy: null,
      DTM: null,

      // Filter topic
      current_topic: -2,
      topics: new Set(),

      //Summirise
      summarise: null,

      //Graphs
      messageGraph: null,
      messageStatusGraph: 200,

      current_graph: {label: 'authors'},
      list_of_graphs: [{label: 'authors'}, {label: 'journals'}, {label: 'countries'}],

      infoAuthorsData: null,
      infoJournalsData: null,
      infoCountryData: null
    }
  }


  componentDidMount() {
    this.getArticles();
    this.getGraphInfo();
    console.log('start');
  }

  // Search

  getArticles = (url, interval = 1000) => {
    fetch(variables.API_URL + '/api/search/',
      {
        headers: {
          'Content-Type': 'application/json;charset=utf-8',
          'Authorization': `Token ${variables.token}`,
        },
      }
    )
      .then(response => {
        console.log(response.status);
        ErrorMessage = response.status
        if (response.ok) {
          return response.json()
        } else {
          throw Error(response.status)
        }
      })
      .then(data => {
        if (ErrorMessage === 202) {
          this.setState({ loading: true, message: data.message, messageStatus: 202 });
          console.log(data.message);
          setTimeout(() => {
            return this.getArticles(url, interval)
          }, interval);
        } else {
            console.log('This work')
            console.log(data.task.query.split('AND ')[1].split(':')[1].replace('/', '-').replace('/', '-'))
            this.setState({
              articles: data.data.search_ncbi,
              loading: false,
              full_query: data.task.full_query,
              translation_stack: data.task.translation_stack,
              short_query: data.task.query,
              count: data.task.count,
              message: "Запрос успешно обработан",
              messageStatus: 200,
              queryStartDate: data.task.query.split('AND ')[1].split(':')[0].replace('/', '-').replace('/', '-'),
              queryEndDate: data.task.query.split('AND ')[1].split(':')[1].replace('/', '-').replace('/', '-').replace('[dp]', ''),
              queryText: data.task.query.split(' AND')[0]

            });
        }
      })
      .catch(error => {
        if (ErrorMessage === 500) {
            this.setState({ articles: [], DetailArticle: null, loading: false, message: 'Ошибка сервера', messageStatus: 500 });
        } else {
            this.setState({ articles: [], DetailArticle: null, loading: false, message: 'Что-то пошло не так', messageStatus: 400 });
        }
      })
  }

  createTask() {
    // Отправляем запрос на сервер для получения статей
    if (this.state.queryText === '') {
        alert('Пожайлуста заполните поле запроса!');
        return;
    }
    fetch(variables.API_URL + '/api/search/',
      {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json;charset=utf-8',
          'Authorization': `Token ${this.state.token}`,
        },
        body: JSON.stringify({
          search_field: this.state.queryText,
          dateStart: this.state.queryStartDate,
          dateStop: this.state.queryEndDate,
          Gender: [...this.state.queryGenders],
          Type: [...this.state.queryTypes],
          Old: [...this.state.queryOlds],
        })
      })
      .then(response => {
        console.log(response.status);
        if (response.ok) {
          return response.json()
        } else {
          ErrorMessage = response.status
          throw Error(response.status)
        }
      })
      .then(data => {
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
        this.getArticles()
      })
      .catch(error => {
        console.log(error)
        if (ErrorMessage === 500) {
            this.setState({ task: null, loading: false, message: 'Ошибка сервера', messageStatus: 500 });
        } else if (ErrorMessage === 403) {
            this.setState({ task: null, loading: false, message: 'Дождитесь окончания предыдушего запроса', messageStatus: 400 });
        } else {
            this.setState({ task: null, loading: false, message: 'Что=то пошло не так', messageStatus: 400 });
        }
      })
  }

  componentDidMount() {
    this.getArticles();
    this.getAnalise();
    this.getGraphInfo();
    console.log('start search');
  }

  onSelectionAnalise = () => {
    const selectedRows = this.gridRef.current.api.getSelectedRows();
    this.setState({ DetailArticle: (selectedRows.length === 1 ? selectedRows[0] : null) })
  }

  changeQueryText = (e) => {
    this.setState({ queryText: e.target.value });
  }

  changeQueryStartDate = (e) => {
    this.setState({ queryStartDate: e.target.value });
  }

  changeQueryEndDate = (e) => {
    this.setState({ queryEndDate: e.target.value });
  }

  changeQueryGenders(gender) {
    if (this.state.queryGenders.has(gender)) {
      this.state.queryGenders.delete(gender)
    } else {
      this.state.queryGenders.add(gender)
    }
    this.setState({ updateOr: !this.state.updateOr })

  }

  changeQueryTypes(type) {
    if (this.state.queryTypes.has(type)) {
      this.state.queryTypes.delete(type)
    } else {
      this.state.queryTypes.add(type)
    }
    this.setState({ updateOr: !this.state.updateOr })
  }

  changeQueryOlds(old) {
    if (this.state.queryOlds.has(old)) {
      this.state.queryOlds.delete(old)
    } else {
      this.state.queryOlds.add(old)
    }
    this.setState({ updateOr: !this.state.updateOr })
  }

  changeQueryOldsMany(olds) {
    if (!this.state.useAll) {
      for (let old of olds) {
        this.state.queryOlds.delete(old)
      }
    } else {
      for (let old of olds) {
        this.state.queryOlds.add(old)
      }
    }
    this.setState({ useAll: !this.state.useAll })
  }

  RoundPersent(number) {
    return number.toFixed(2);
  }

  startAnalise() {
    let analise_data = [];
    this.gridRef.current.api.forEachNodeAfterFilter((rowNode) => analise_data.push(rowNode.data.uid));
    fetch(variables.API_URL + '/api/analise/', {
      method: 'POST',
      headers: {
        'Accept': 'application/json',
        'Content-Type': 'application/json;charset=utf-8',
        'Authorization': `Token ${variables.token}`,
      },
      body: JSON.stringify({
        articles: analise_data
      })
    })
      .then(response => {
        console.log(response.status);
        if (response.ok) {
          return response.json()
        } else {
          ErrorMessage = response.status
          throw Error(response.status)
        }
      })
      .then((data) => {
        this.setState({ messageAnalise: "Ваш запрос в очереди. Пожайлуста дождитесь результата", messageStatusAnalise: 201 })
        this.getAnalise();
      })
      .catch((error) => {
        if (ErrorMessage === 500) {
            this.setState({ data: [], dataInfo: [], DetailArticle: null, loading: false, messageAnalise: 'Ошибка сервера', messageStatusAnalise: 500 });
        } else if (ErrorMessage === 403) {
            this.setState({ data: [], dataInfo: [], DetailArticle: null, loading: false, messageAnalise: 'Дождитесь окончания предыдушего запроса', messageStatusAnalise: 403 });
        } else {
            this.setState({ data: [], dataInfo: [], DetailArticle: null, loading: false, messageAnalise: 'Что=то пошло не так', messageStatusAnalise: 400 });
        }
      })
  }

  // Analise

  onSelectionChanged = (gridApi) => {
    const selectedRows = this.gridAnaliseRef.current.api.getSelectedRows();
    this.setState({ DetailArticle: (selectedRows.length === 1 ? selectedRows[0] : null) })
  }

  getAnalise = (url, interval = 1000) => {
    fetch(variables.API_URL + "/api/analise/",
      {
        headers: {
          'Content-Type': 'application/json;charset=utf-8',
          'Authorization': `Token ${variables.token}`,
        },
      }
    )
      .then((res) => {
        if (res.ok) {
          ErrorMessage = res.status
          return res.json()
        } else {
          throw new Error(res)
        }
      })
      .then((data) => {
        console.log(data)
        if (ErrorMessage === 202) {
          this.setState({ loading: true, messageAnalise: data.message, messageStatusAnalise: 202 });
          setTimeout(() => {
            return this.getAnalise(url, interval)
          }, interval);
        } else {
          if (data.data.clust_graph !== null) {
            // delete data.data.clust_graph.layout.width;
            data.data.clust_graph.layout.width = 800;
          }
          if (data.data.heapmap !== null)
            if (data.data.heapmap.length !== 0) {
                delete data.data.heapmap.layout.width;
            }
          if (data.data.heirarchy !== null)
            if (data.data.heirarchy.length !== 0) {
                delete data.data.heirarchy.layout.width;
            }
          if (data.data.DTM !== null)
            if (data.data.DTM.length !== 0) {
                delete data.data.DTM.layout.width;
            }
          var topics = new Set()
          if (data.data.tematic_analise.length !== 0) {
              for (let record of data.data.tematic_analise) {
                if (!topics.has(record.topic)) {
                  topics.add(record.topic)
                }
              }
          }
          console.log(topics)
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
            topics: [...topics],
          });
          this.getGraphInfo()
        }
      })
      .catch((error) => {
        console.log(error);
        if (ErrorMessage === 500) {
            this.setState({ analise_articles: [], clust_graph: null, heapmap: null, heirarchy: null, DetailArticle: null, loading: false, messageAnalise: 'Ошибка сервера', messageStatusAnalise: 500 });
        } else if (ErrorMessage === 403) {
            this.setState({ data: [], dataInfo: [], DetailArticle: null, loading: false, messageAnalise: 'Дождитесь окончания предыдушего запроса', messageStatusAnalise: 403 });
        } else {
            this.setState({ analise_articles: [], clust_graph: null, heapmap: null, heirarchy: null, DetailArticle: null, loading: false, messageAnalise: 'Что-то пошло не так', messageStatusAnalise: 400 });
        }
      });
  }

  externalFilterChanged = (newValue) => {
    this.setState({ current_topic: newValue.x, summarise: null });
    this.gridAnaliseRef.current.api.onFilterChanged();
  }

  isExternalFilterPresent = () => {
    // if ageType is not everyone, then we are filtering
    return this.state.current_topic !== -2;
  }

  doesExternalFilterPass = (node) => {
    if (node.data) {
      if (node.data.topic === this.state.current_topic) {
        return true
      }
    }
    return false;
  }

  autoGroupColumnDef = () => {
    return {
      minWidth: 200,
    }
  }

  getGraphData = () => {
    var current_topic = this.state.current_topic.toString();
    if (current_topic === '-2') {
      return this.state.clust_graph.data
    }

    var data = []
    for (let topic of this.state.clust_graph.data) {
      var topic_id = topic.name.split('_')[0];
      if (current_topic === topic_id) {
        data.push(topic);
      }
    }
    return data
  }

  getSummarise = (task_id, interval = 1000) => {
    fetch(variables.API_URL + `/api/summarise?task_id=${task_id}`,
      {
        headers: {
          'Content-Type': 'application/json;charset=utf-8',
          'Authorization': `Token ${variables.token}`,
        },
      }
    )
      .then((res) => {
        if (res.status == 202) {
          this.setState({ loading: true, messageStatus: 202, message: "Отправлено на суммаризацию..." })
          setTimeout(() => {
            return this.getSummarise(task_id, interval)
          }, interval);
        }
        if (res.status == 200) {
          return res.json()
        } else {
          ErrorMessage = res.status
          throw Error(res.statusText)
        }
      })
      .then((data) => {
        this.setState({
          summarise: data.data,
           loading: false,
          message: 'Суммаризация прошла успешно',
          messageStatus: 200
        });
      })
      .catch((err) => {
        console.log(err);
        this.setState({ summarise: null, message: 'Произошла ошибка при суммаризации', messageStatus: 500 });
      });
  }

  createSummariseQuery() {
    if (this.state.current_topic === -2) {
      var data = this.state.analise_articles
    } else {
      var data = []
      for (let article of this.state.analise_articles) {
        if (this.state.current_topic === article.topic) {
          data.push(article);
        }
      }
    }
    fetch(variables.API_URL + '/api/summarise', {
      method: 'POST',
      headers: {
        'Accept': 'application/json',
        'Content-Type': 'application/json;charset=utf-8',
        'Authorization': `Token ${variables.token}`,
      },
      body: JSON.stringify({
        articles: data
      })
    })
      .then(response => {
        console.log(response.status);
        if (response.ok) {
          return response.json()
        } else {
          ErrorMessage = response.status
          throw Error(response.status)
        }
      })
      .then((result) => {
        var task_id = result.data;
        this.setState({ message: 'Отправлено на суммаризацию пожайлуста дождитесь ответа', messageStatus: 201 })
        this.getSummarise(task_id);
      })
      .catch((error) => {
        this.setState({ message: 'Ошибка при суммаризации', messageStatus: 500 })
      })
  }

  // Graphs authors, countries, jornals
  getGraphInfo() {
    fetch(variables.API_URL + '/api/graphs/', {
      headers: {
        'Accept': 'application/json',
        'Content-Type': 'application/json;charset=utf-8',
        'Authorization': `Token ${variables.token}`,
      }
    })
      .then(response => {
        console.log(response.status);
        if (response.ok) {
          return response.json()
        } else {
          ErrorMessage = response.status
          throw Error(response.status)
        }
      })
      .then((data) => {
        this.setState({ infoAuthorsData: data.data.info_graph, infoJournalsData: data.data.info_graph_journals, infoCountryData: data.data.info_graph_countries})
      })
      .catch((error) => {
        console.log(error)
      })
  }

  // разметка

  getMarkUp = (task_id, interval = 1000) => {
    fetch(variables.API_URL + `/api/markup?task_id=${task_id}`,
      {
        headers: {
          'Content-Type': 'application/json;charset=utf-8',
          'Authorization': `Token ${variables.token}`,
        },
      }
    )
      .then((res) => {
        if (res.status == 202) {
          setTimeout(() => {
            return this.getMarkUp(task_id, interval)
          }, interval);
        } else if (res.status == 200) {
          return res.json()
        } else {
          throw Error(res.statusText)
        }
      })
      .then((data) => {
        try {
          this.setState({
            DetailArticle: data.data,
            message: 'Разметка прошла успешно',
            messageStatus: 200,
            loading: false,
          });
        } catch {
          console.log('access')
        }
      })
      .catch((err) => {
        console.log(err);
        this.setState({ message: 'Произошла ошибка при разметке', loading: false, messageStatus: 500 });
      });
  }

  markUpArticle(DetailArticle) {
    this.setState({ loading: true })
    fetch(variables.API_URL + '/api/markup', {
      method: 'POST',
      headers: {
        'Accept': 'application/json',
        'Content-Type': 'application/json;charset=utf-8',
        'Authorization': `Token ${variables.token}`,
      },
      body: JSON.stringify({
        article: DetailArticle
      })
    })
      .then((res) => {
        console.log(res.status)
        if (res.ok) {
          return res.json()
        } else {
          throw Error(res.statusText)
        }
      })
      .then((result) => {
        var task_id = result.data;
        this.setState({ message: 'Отправлено на суммаризацию пожайлуста дождитесь ответа', messageStatus: 201 })
        this.getMarkUp(task_id);
      })
      .catch((err) => {
        console.log(err);
        this.setState({
          message: 'ошибка при разметке',
          messageStatus: 500,
          loading: false,
        });
      });
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

      analise_articles,
      analise_info,
      clust_graph,
      heapmap,
      heirarchy,
      DTM,
      current_topic,
      topics,
      summarise,
      allow_page,

      current_graph,
      list_of_graphs,
      infoAuthorsData,
      infoCountryData,
      infoJournalsData
    } = this.state;

    if (!token) {
      return <Navigate push to="/login" />
    } else if (allow_page === 1) {
      return <Navigate push to="/ddi_review" />
    } else {
      return (
        <>
          <header>
            <nav class="bg-white px-4 lg:px-6 py-2.5">
              <div class="flex flex-wrap justify-between items-center">
                <div class="flex justify-start items-center">
                  <a href="" class="flex mr-4">
                    <img src="https://flowbite.s3.amazonaws.com/logo.svg" class="mr-3 h-8" alt="FlowBite Logo" />
                    <span class="self-center text-2xl font-semibold whitespace-nowrap">EBM DаtaMed</span>
                  </a>
                  {allow_page === 2?
                      <ul class="flex font-medium flex-row space-x-8 ml-10">
                        <Link to="/tematic_review">
                          <li>
                            <a href="#" class="block py-2 pl-3 pr-4 text-gray-900 bg-blue-700 rounded md:bg-transparent md:text-blue-700 md:p-0" aria-current="page">Тематический анализ</a>
                          </li>
                        </Link>
                        <Link to="/ddi_review">
                          <li>
                            <a href="#" class="block py-2 pl-3 pr-4 text-gray-900 rounded hover:bg-gray-100 md:hover:bg-transparent md:hover:text-blue-700 md:p-0">Факты для EBM</a>
                          </li>
                        </Link>
                      </ul>
                  :null}

                </div>
                <div class="flex items-center lg:order-2">
                  <div class="flex-shrink-0 dropdown">
                    <a href="#" class="d-block link-body-emphasis text-decoration-none dropdown-toggle" data-bs-toggle="dropdown" aria-expanded="false">
                      <img src="https://github.com/mdo.png" alt="mdo" width="32" height="32" class="rounded-circle" />
                    </a>
                    <ul class="dropdown-menu text-small shadow">
                      <li><a class="dropdown-item" href="#">New project...</a></li>
                      <li><a class="dropdown-item" href="#">Settings</a></li>
                      <li><a class="dropdown-item" href="#">Profile</a></li>
                      <li><a class="dropdown-item" href="#">Sign out</a></li>
                    </ul>
                  </div>
                </div>
              </div>
            </nav>
            <nav class="bg-white border border-gray-200 py-2 px-6">
              <div class="w-full">
                <div class="flex justify-between items-center">
                  <button id="toggleSidebar" aria-expanded="true" aria-controls="sidebar" class="hidden order-first p-2 mr-3 text-gray-600 rounded cursor-pointer lg:inline hover:text-gray-900 hover:bg-gray-100 dark:text-gray-400 dark:hover:text-white dark:hover:bg-gray-700" data-bs-toggle="collapse" data-bs-target="#sidebar" aria-label="Toggle navigation">
                    <svg class="w-6 h-6" fill="currentColor" viewBox="0 0 20 20" xmlns="http://www.w3.org/2000/svg"><path fill-rule="evenodd" d="M3 5a1 1 0 011-1h12a1 1 0 110 2H4a1 1 0 01-1-1zM3 10a1 1 0 011-1h6a1 1 0 110 2H4a1 1 0 01-1-1zM3 15a1 1 0 011-1h12a1 1 0 110 2H4a1 1 0 01-1-1z" clip-rule="evenodd"></path></svg>
                  </button>
                  <label for="topbar-search" class="sr-only">Поисковый запрос</label>
                  <div>
                    <label for="search" class="mb-2 text-sm font-medium text-gray-900 sr-only dark:text-white">Поисковый запрос</label>
                    <div class="relative mt-1 grow lg:w-96">
                      <div class="absolute inset-y-0 left-0 flex items-center pl-3 pointer-events-none">
                        <svg class="w-4 h-4 text-gray-500 dark:text-gray-400" aria-hidden="true" xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 20 20">
                          <path stroke="currentColor" stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="m19 19-4-4m0-7A7 7 0 1 1 1 8a7 7 0 0 1 14 0Z" />
                        </svg>
                      </div>
                      <input
                        class="py-3 bg-gray-50 border border-gray-300 text-gray-900 sm:text-sm rounded-lg focus:ring-primary-500 focus:border-primary-500 block w-full pl-10 p-2.5"
                        id="search"
                        type="text"
                        name="search_field"
                        placeholder={short_query? short_query.split(' AND')[0]: 'covid-19'}
                        value={queryText}
                        onChange={this.changeQueryText}
                        aria-label="Search" />
                      <button type="submit" value="Найти" onClick={() => this.createTask()} class="text-white absolute right-2.5 bottom-2.5 bg-blue-700 hover:bg-blue-800 focus:ring-4 focus:outline-none focus:ring-blue-300 font-medium rounded-lg text-sm px-4 py-2">Найти</button>
                    </div>
                  </div>
                  <div class="ml-5 grow items-center justify-between hidden w-full md:flex md:w-auto md:order-1" id="navbar-sticky">
                    <ul class="nav nav-pills" id="myTab" role="tablist">
                      <li class="nav-item mr-2" role="presentation">
                        <button class="nav-link inline-block px-4 py-2 text-white bg-blue-600 rounded-lg active" id="home-tab" data-bs-toggle="tab" data-bs-target="#home" type="button" role="tab" aria-controls="home" aria-selected="false">Результаты поиска</button>
                      </li>
                      <li class="nav-item mr-2" role="presentation">
                        <button class="nav-link inline-block px-4 py-2 rounded-lg hover:text-gray-900 hover:bg-gray-100" id="profile-tab" data-bs-toggle="tab" data-bs-target="#profile" type="button" role="tab" aria-controls="profile" aria-selected="true">Тематическое описание коллекции</button>
                      </li>
                      <li class="nav-item mr-2" role="presentation">
                        <button class="nav-link inline-block px-4 py-2 rounded-lg hover:text-gray-900 hover:bg-gray-100" id="contact-tab" data-bs-toggle="tab" data-bs-target="#contact" type="button" role="tab" aria-controls="contact" aria-selected="false" >Схема</button>
                      </li>
                    </ul>
                    <button id="toggleSidebar" aria-expanded="true" aria-controls="sidebar2" class="order-last hidden p-2 text-gray-600 rounded cursor-pointer lg:inline hover:text-gray-900 hover:bg-gray-100 dark:text-gray-400 dark:hover:text-white dark:hover:bg-gray-700" data-bs-toggle="collapse" data-bs-target="#sidebar2" aria-label="Toggle navigation">
                      <svg class="w-6 h-6 rotate-180" fill="currentColor" viewBox="0 0 20 20" xmlns="http://www.w3.org/2000/svg"><path fill-rule="evenodd" d="M3 5a1 1 0 011-1h12a1 1 0 110 2H4a1 1 0 01-1-1zM3 10a1 1 0 011-1h6a1 1 0 110 2H4a1 1 0 01-1-1zM3 15a1 1 0 011-1h12a1 1 0 110 2H4a1 1 0 01-1-1z" clip-rule="evenodd"></path></svg>
                    </button>
                  </div>
                </div>
              </div>
            </nav>
          </header >
          <main>
            <div>
              <div className="container-fluid h-screen">
                <div className="row align-items-stretch b-height">
                  <aside id="sidebar" className="h-screen col-md-2 my-3 bg-white collapse show width border rounded-3 overflow-auto g-0">
                    <div className="accordion accordion-flush" id="accordionFlushExample">
                      <div className="accordion-item">
                        <h2 className="accordion-header" id="flush-headingOne">
                          <button className="accordion-button collapsed" type="button" data-bs-toggle="collapse" data-bs-target="#flush-collapseOne" aria-expanded="false" aria-controls="flush-collapseOne">
                            Дата публикации
                          </button>
                        </h2>
                        <div id="flush-collapseOne" className="collapse show multi-collapse" aria-labelledby="flush-headingOne" data-bs-target="#accordionFlushExample">
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
                          <button class="accordion-button collapsed" data-target='#flush-collapseThree' type="button" data-bs-toggle="collapse" data-bs-target="#flush-collapseFour" aria-expanded="false" aria-controls="flush-collapseThree">
                            Тип статьи
                          </button>
                        </h2>
                        <div id="flush-collapseFour" className="collapse show multi-collapse" aria-labelledby="flush-headingFour" data-bs-target="#accordionFlushExample">
                          <div className="accordion-body">
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxClinicalTrial"
                                name="CheckBoxClinicalTrial"
                                checked={this.state.queryTypes.has('clinical trial')}
                                onChange={() => this.changeQueryTypes('clinical trial')}
                              />
                              <label className="form-check-label" for="CheckboxClinicalTrial">Clinical Trial</label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxMetaAnalysys"
                                name="CheckboxMetaAnalysys"
                                checked={this.state.queryTypes.has('meta-analysis')}
                                onChange={() => this.changeQueryTypes('meta-analysis')}
                              />
                              <label className="form-check-label" for="CheckboxMetaAnalysys">Meta Analysys</label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxRandomizedControlledTrial"
                                name="CheckboxRandomizedControlledTrial"
                                checked={this.state.queryTypes.has('randomized controlled trial')}
                                onChange={() => this.changeQueryTypes('randomized controlled trial')}
                              />
                              <label className="form-check-label" for="CheckboxRandomizedControlledTrial">Randomized Controlled Trial</label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxReview"
                                name="CheckboxReview"
                                checked={this.state.queryTypes.has('review')}
                                onChange={() => this.changeQueryTypes('review')}
                              />
                              <label className="form-check-label" for="CheckboxReview">Review</label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxSystematicReview"
                                name="CheckboxSystematicReview"
                                checked={this.state.queryTypes.has('systematic review')}
                                onChange={() => this.changeQueryTypes('systematic review')}
                              />
                              <label className="form-check-label" for="CheckboxSystematicReview">Systematic Review</label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxJournalArticle"
                                name="CheckboxJournalArticle"
                                checked={this.state.queryTypes.has('journal article')}
                                onChange={() => this.changeQueryTypes('journal article')}
                              />
                              <label className="form-check-label" for="CheckboxJournalArticle">Journal Article</label>
                            </div>

                          </div>
                        </div>
                      </div>
                      <div className="accordion-item">
                        <h2 className="accordion-header" id="flush-headingThree">
                          <button className="accordion-button collapsed" data-target='#flush-collapseThree' type="button" data-bs-toggle="collapse" data-bs-target="#flush-collapseThree" aria-expanded="false" aria-controls="flush-collapseThree">
                            Возраст пациента
                          </button>
                        </h2>
                        <div id="flush-collapseThree" className="collapse show multi-collapse" aria-labelledby="flush-headingThree" data-bs-target="#accordionFlushExample">
                          <div className="accordion-body">
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxChild18"
                                checked={this.state.queryOlds.has('infant[mh]') & this.state.queryOlds.has('child[mh]') & this.state.queryOlds.has('adolescent[mh]')}
                                onChange={() => this.changeQueryOldsMany(['infant[mh]', 'child[mh]', 'adolescent[mh]'])}
                              />
                              <label className="form-check-label" for="CheckboxChild18">Child: birth-18 years</label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxNewborn"
                                checked={this.state.queryOlds.has('infant, newborn[mh]')}
                                onChange={() => this.changeQueryOlds('infant, newborn[mh]')}
                              />
                              <label className="form-check-label" for="CheckboxNewborn">Newborn: birth-1 months</label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxInfant023"
                                checked={this.state.queryOlds.has('infant[mh]')}
                                onChange={() => this.changeQueryOlds('infant[mh]')}
                              />
                              <label className="form-check-label" for="CheckboxInfant023">Infant: birth-23 months</label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxInfant123"
                                checked={this.state.queryOlds.has('infant[mh:noexp]')}
                                onChange={() => this.changeQueryOlds('infant[mh:noexp]')}
                              />
                              <label className="form-check-label" for="CheckboxInfant123">Infant: 1-23 months</label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxPreschool"
                                checked={this.state.queryOlds.has('child, preschool[mh]')}
                                onChange={() => this.changeQueryOlds('child, preschool[mh]')}
                              />
                              <label className="form-check-label" for="CheckboxPreschool">Preschool Child: 2-5 years</label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxChild612"
                                checked={this.state.queryOlds.has('child[mh:noexp]')}
                                onChange={() => this.changeQueryOlds('child[mh:noexp]')}
                              />
                              <label className="form-check-label" for="CheckboxChild612">Child: 6-12 years</label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxAdolescent"
                                checked={this.state.queryOlds.has('adolescent[mh]')}
                                onChange={() => this.changeQueryOlds('adolescent[mh]')}
                              />
                              <label className="form-check-label" for="CheckboxAdolescent">Adolescent: 13-18 years</label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxAdult19"
                                checked={this.state.queryOlds.has('adult[mh]')}
                                onChange={() => this.changeQueryOlds('adult[mh]')}
                              />
                              <label className="form-check-label" for="CheckboxAdult19">Adult: 19+ years</label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxYoungAdult"
                                checked={this.state.queryOlds.has('young adult[mh]')}
                                onChange={() => this.changeQueryOlds('young adult[mh]')}
                              />
                              <label className="form-check-label" for="CheckboxYoungAdult">Young Adult: 19-24 years</label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxAdult1944"
                                checked={this.state.queryOlds.has('adult[mh:noexp]')}
                                onChange={() => this.changeQueryOlds('adult[mh:noexp]')}
                              />
                              <label className="form-check-label" for="CheckboxAdult1944">Adult: 19-44 years</label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxMiddle45"
                                checked={this.state.queryOlds.has('middle aged[mh]') & this.state.queryOlds.has('aged[mh]')}
                                onChange={() => this.changeQueryOldsMany(['middle aged[mh]', 'aged[mh]'])}
                              />
                              <label className="form-check-label" for="CheckboxMiddle45">Middle Aged: + Aged 45+ years</label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxMiddle4565"
                                checked={this.state.queryOlds.has('middle aged[mh]')}
                                onChange={() => this.changeQueryOlds('middle aged[mh]')}
                              />
                              <label className="form-check-label" for="CheckboxMiddle4565">Middle Aged: 45-65 years</label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxAged"
                                checked={this.state.queryOlds.has('aged[mh]')}
                                onChange={() => this.changeQueryOlds('aged[mh]')}
                              />
                              <label className="form-check-label" for="CheckboxAged">Aged: 65+ years</label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="Checkbox80"
                                checked={this.state.queryOlds.has('aged, 80 and over[mh]')}
                                onChange={() => this.changeQueryOlds('aged, 80 and over[mh]')}
                              />
                              <label className="form-check-label" for="Checkbox80">80 and over: 80+ years</label>
                            </div>
                          </div>
                        </div>
                      </div>
                      <div className="accordion-item">
                        <h2 className="accordion-header" id="flush-headingTwo">
                          <button className="accordion-button collapsed" type="button" data-bs-toggle="collapse" data-bs-target="#flush-collapseTwo" aria-expanded="false" aria-controls="flush-collapseTwo">
                            Пол
                          </button>
                        </h2>
                        <div id="flush-collapseTwo" className="collapse show multi-collapse" aria-labelledby="flush-headingTwo" data-bs-target="#accordionFlushExample">
                          <div className="accordion-body">
                            <div className="form-check form-check-inline">
                              <input
                                type="checkbox"
                                id="CheckboxMan"
                                name="CheckboxMan"
                                className="form-check-input"
                                checked={this.state.queryGenders.has('man')}
                                onChange={() => this.changeQueryGenders('man')}
                              />
                              <label className="form-check-label" for="CheckboxMan">Мужчина</label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                name="CheckboxWomen"
                                id="CheckboxWomen"
                                checked={this.state.queryGenders.has('woman')}
                                onChange={() => this.changeQueryGenders('woman')}
                              />
                              <label className="form-check-label" for="CheckboxWomen">Женщина</label>
                            </div>
                          </div>
                        </div>
                      </div>
                    </div>
                  </aside>
                  <section class="col p-3 m-3 border rounded-3 bg-white h-screen overflow-auto">
                    <div class="bd-example">
                      <div class="tab-content" id="myTabContent">
                        <div class="tab-pane fade active show" id="home" role="tabpanel" aria-labelledby="home-tab">
                          <div class="container-fluid g-0">
                            <div class="accordion accordion-flush" id="accordion">
                              <div class="accordion-item">
                                <h2 class="accordion-header" id="">
                                  {message?
                                      messageStatus > 299?
                                        <p class="pb-2 mb-3 border-bottom" style={{ color: 'red' }}>{message}. </p>
                                        : messageStatus === 200 ?
                                        <p class="pb-2 mb-3 border-bottom" style={{ color: 'green' }}>{message}.</p>
                                        :
                                        <p class="pb-2 mb-3 border-bottom" style={{ color: 'black' }}>{message}.</p>
                                      :null
                                  }
                                  <button class="accordion-button collapsed" type="button" data-bs-toggle="collapse" data-bs-target="#flush-collapseEleven" aria-expanded="false" aria-controls="flush-collapseSeven">
                                    По запросу найдено {count} источников.
                                  </button>
                                </h2>
                                <div id="flush-collapseEleven" class="collapse multi-collapse" aria-labelledby="flush-headingEleven" data-bs-target="#accordionFlushExample">
                                  <div class="accordion-body">
                                    <p class="pb-2 mb-3 border-bottom"> Запрос {short_query} .</p>
                                    <p class="pb-2 mb-3 border-bottom"> Запрос автоматически расширен до следующего вида - {full_query}.</p>
                                    <p class="pb-2 mb-3 border-bottom"> Служебная информация для анализа : {translation_stack}.</p>
                                  </div>
                                </div>
                              </div>
                            </div>
                            <input className="text-white right-2.5 my-4 bottom-2.5 bg-blue-700 hover:bg-blue-800 focus:ring-4 focus:outline-none focus:ring-blue-300 font-medium rounded-lg text-sm px-4 py-2" type="submit" value="Начать обработку" onClick={() => this.startAnalise()} />

                            <div className="ag-theme-alpine ag-theme-acmecorp" style={{ height: 700 }}>
                              <AgGridReact
                                ref={this.gridRef}
                                rowData={articles}
                                columnDefs={articlesInfo}
                                pagination={true}
                                onSelectionChanged={this.onSelectionAnalise}
                                autoGroupColumnDef={this.autoGroupColumnDef}
                                onChange={this.externalFilterChanged}
                                rowSelection={'single'}
                                sideBar={{
                                  toolPanels: [
                                    {
                                      id: 'columns',
                                      labelDefault: 'Columns',
                                      labelKey: 'columns',
                                      iconKey: 'columns',
                                      toolPanel: 'agColumnsToolPanel',
                                      minWidth: 225,
                                      width: 225,
                                      maxWidth: 225,
                                    },
                                    {
                                      id: 'filters',
                                      labelDefault: 'Filters',
                                      labelKey: 'filters',
                                      iconKey: 'filter',
                                      toolPanel: 'agFiltersToolPanel',
                                      minWidth: 180,
                                      maxWidth: 400,
                                      width: 250,
                                    },
                                  ],
                                  position: 'left'
                                }}
                              >
                              </AgGridReact>
                            </div>
                          </div>
                        </div>
                        <div class="tab-pane fade" id="profile" role="tabpanel" aria-labelledby="profile-tab">
                          <div>
                            <h2 class="accordion-header" id="">
                                  {messageAnalise?
                                      messageStatusAnalise > 299?
                                        <p class="pb-2 mb-3 border-bottom" style={{ color: 'red' }}>{messageAnalise}. </p>
                                        : messageStatusAnalise === 200 ?
                                        <p class="pb-2 mb-3 border-bottom" style={{ color: 'green' }}>{messageAnalise}.</p>
                                        :
                                        <p class="pb-2 mb-3 border-bottom" style={{ color: 'black' }}>{messageAnalise}.</p>
                                      :null
                                  }
                            </h2>
                            <p>Выбрана тема: {current_topic === -2 ? "Выбраны все" : `№ ${current_topic} из ${topics.length}`}</p>
                            <Slider
                              axis="x"
                              x={current_topic}
                              xmax={topics.length - 2}
                              xmin={-2}
                              onChange={this.externalFilterChanged}
                            />
                          </div>
                          <div class="accordion accordion-flush" id="accordion">
                              <div className="accordion-item">
                                <h2 className="accordion-header" id="flush-headingTwo">
                                  <button className="accordion-button collapsed" type="button" data-bs-toggle="collapse" data-bs-target="#flush-collapseSix" aria-expanded="false" aria-controls="flush-collapseTwo">
                                    Таблица
                                  </button>
                                </h2>
                                <div id="flush-collapseSix" className="collapse show multi-collapse" aria-labelledby="flush-headingTwo" data-bs-target="#accordionFlushExample">
                                  <div className="accordion-body">
                                    <div className="ag-theme-alpine ag-theme-acmecorp" style={{ height: 700 }}>
                                      <AgGridReact
                                        ref={this.gridAnaliseRef}
                                        rowData={analise_articles}
                                        columnDefs={analise_info}
                                        pagination={true}
                                        rowSelection={'single'}
                                        onSelectionChanged={this.onSelectionChanged}
                                        animateRows={true}
                                        isExternalFilterPresent={this.isExternalFilterPresent}
                                        doesExternalFilterPass={this.doesExternalFilterPass}
                                        sideBar={{
                                          toolPanels: [
                                            {
                                              id: 'columns',
                                              labelDefault: 'Columns',
                                              labelKey: 'columns',
                                              iconKey: 'columns',
                                              toolPanel: 'agColumnsToolPanel',
                                              minWidth: 225,
                                              width: 225,
                                              maxWidth: 225,
                                            },
                                            {
                                              id: 'filters',
                                              labelDefault: 'Filters',
                                              labelKey: 'filters',
                                              iconKey: 'filter',
                                              toolPanel: 'agFiltersToolPanel',
                                              minWidth: 180,
                                              maxWidth: 400,
                                              width: 250,
                                            },
                                          ],
                                          position: 'left',
                                        }}
                                      >
                                      </AgGridReact>
                                    </div>
                                  </div>
                                </div>
                              </div>
                          </div>
                          <div class="accordion accordion-flush" id="accordion">
                            <div className="accordion-item">
                              <h2 className="accordion-header" id="flush-headingTwo">
                                <button className="accordion-button collapsed" type="button" data-bs-toggle="collapse" data-bs-target="#flush-collapseSeven" aria-expanded="false" aria-controls="flush-collapseTwo">
                                  Класстеризация
                                </button>
                              </h2>
                              <div id="flush-collapseSeven" className="collapse show multi-collapse" aria-labelledby="flush-headingTwo" data-bs-target="#accordionFlushExample">
                                <div className="accordion-body">
                                  <div>
                                    {clust_graph ?
                                      <Plot
                                        data={this.getGraphData()}
                                        layout={clust_graph.layout}
                                      />
                                      : null
                                    }
                                  </div>
                                </div>
                              </div>
                            </div>

                            <div className="accordion-item">
                              <h2 className="accordion-header" id="flush-headingTwo">
                                <button className="accordion-button collapsed" type="button" data-bs-toggle="collapse" data-bs-target="#flush-collapseEight" aria-expanded="false" aria-controls="flush-collapseTwo">
                                  Таблица близости
                                </button>
                              </h2>
                              <div id="flush-collapseEight" className="collapse show multi-collapse" aria-labelledby="flush-headingTwo" data-bs-target="#accordionFlushExample">
                                <div className="accordion-body">
                                  <div>
                                    {heapmap ?
                                      <Plot
                                        data={heapmap.data}
                                        layout={heapmap.layout}
                                      />
                                      : null
                                    }
                                  </div>
                                </div>
                              </div>
                            </div>
                            <div className="accordion-item">
                              <h2 className="accordion-header" id="flush-headingTwo">
                                <button className="accordion-button collapsed" type="button" data-bs-toggle="collapse" data-bs-target="#flush-collapseNine" aria-expanded="false" aria-controls="flush-collapseTwo">
                                  Иерархическое дерево
                                </button>
                              </h2>
                              <div id="flush-collapseNine" className="collapse show multi-collapse" aria-labelledby="flush-headingTwo" data-bs-target="#accordionFlushExample">
                                <div className="accordion-body">
                                  <div>
                                    {heirarchy ?
                                      <Plot
                                        data={heirarchy.data}
                                        layout={heirarchy.layout}
                                      />
                                      : null
                                    }
                                  </div>
                                </div>
                              </div>
                            </div>
                            <div className="accordion-item">
                              <h2 className="accordion-header" id="flush-headingTwo">
                                <button className="accordion-button collapsed" type="button" data-bs-toggle="collapse" data-bs-target="#flush-collapseTen" aria-expanded="false" aria-controls="flush-collapseTwo">
                                  Распределение по датам
                                </button>
                              </h2>
                              <div id="flush-collapseTen" className="collapse show multi-collapse" aria-labelledby="flush-headingTwo" data-bs-target="#accordionFlushExample">
                                <div className="accordion-body">
                                  <div>
                                    {DTM ?
                                      <Plot
                                        data={DTM.data}
                                        layout={DTM.layout}
                                      />
                                      : null
                                    }
                                  </div>
                                </div>
                              </div>
                            </div>
                            <div className="accordion-item">
                              <h2 className="accordion-header" id="flush-headingTwo">
                                <button className="accordion-button collapsed" type="button" data-bs-toggle="collapse" data-bs-target="#flush-collapseEleven" aria-expanded="false" aria-controls="flush-collapseTwo">
                                  Суммаризация
                                </button>
                              </h2>
                              <div id="flush-collapseEleven" className="collapse show multi-collapse" aria-labelledby="flush-headingTwo" data-bs-target="#accordionFlushExample">
                                <div className="accordion-body">
                                  <div>
                                    {summarise ?
                                      <>
                                        <p>Summarise</p>
                                        <p>{summarise}</p>
                                      </>
                                      : <input className="btn btn-primary" type="submit" value="Суммаризовать" onClick={() => this.createSummariseQuery()} />}
                                  </div>
                                </div>
                              </div>
                            </div>
                          </div>
                        </div>
                        <div class="tab-pane fade" id="contact" role="tabpanel" aria-labelledby="contact-tab">
                          <div class="container-fluid g-0">
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
                                onChange={(x) => this.setState({current_graph: x})}
                            />
                            <div id="mynetwork" style={{ width: "100%", height: "600px" }}>
                            {current_graph.label === 'authors'?
                                infoAuthorsData ?
                                    <VOSviewerOnline data={infoAuthorsData} />
                                : null
                            :null
                            }
                            {
                            current_graph.label === 'journals'?
                                infoJournalsData ?
                                    <VOSviewerOnline data={infoJournalsData} />
                                : null
                            :null
                            }
                            {current_graph.label === 'countries'?
                                infoCountryData ?
                                    <VOSviewerOnline data={infoCountryData} />
                                : null
                            :null
                            }
                            </div>
                          </div>
                        </div>
                      </div>
                    </div>
                  </section>

                  <aside id="sidebar2" class="col-md-4 h-screen collapse show width col p-3 my-3 border rounded-3 overflow-auto bg-white">
                    <h3 class="pb-2 mb-3 border-bottom">Подробное описание</h3>
                    <nav class="small" id="toc">
                      {DetailArticle ?
                        <div class="card mb-3">
                          <div class="card-body">
                            <a href={DetailArticle.url} class="card-title link-primary text-decoration-none h5"> {DetailArticle.titl} </a>
                            <p class="card-text">---------------------------------- </p>
                            <p class="card-text">Авторы :  {DetailArticle.auth} </p>
                            <p class="card-text">---------------------------------- </p>
                            <p class="card-text">Аннотация :  </p>
                            <p class="card-text" dangerouslySetInnerHTML={{ __html: markup_text(DetailArticle.tiab, DetailArticle.annotations) }} />
                            <p class="card-text">---------------------------------- </p>
                            <p class="card-text"><small class="text-success">Дата публикации : {DetailArticle.pdat} </small></p>
                            <p class="card-text"><small class="text-success">Издание : {DetailArticle.jour}</small></p>
                            <p class="card-text"><small class="text-success">Вид публикации : {DetailArticle.pt}</small></p>
                            <p class="card-text"><small class="text-success">Страна : {DetailArticle.pl} </small></p>
                            <p class="card-text"><small class="text-success">{DetailArticle.mesh} </small></p>
                            {loading ?
                                <p>Loading...</p>
                            : <input className="text-white right-2.5 my-4 bottom-2.5 bg-blue-700 hover:bg-blue-800 focus:ring-4 focus:outline-none focus:ring-blue-300 font-medium rounded-lg text-sm px-4 py-2" type="submit" value="Разметить" onClick={() => this.markUpArticle(DetailArticle)} />
                            }
                          </div>
                        </div>
                        : null}
                    </nav>
                  </aside>

                </div>
              </div>
            </div >
          </main >
        </>
      )
    }
  }
}