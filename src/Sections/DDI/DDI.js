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


import { variables } from '../Variables.js';

import Slider from 'react-input-slider';

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


export class DDIReview extends Component {

  constructor(props) {
    super(props);

    this.gridRef = createRef();
    this.gridAnaliseRef = createRef();
    this.state = {
      token: variables.token,
      loading: false,
      allow_page: variables.allow,

      //queries
      query_list: [],
      message: null,
      messageStatus: 200,
      articles: [],
      articlesInfo: [
        { field: 'text', filter: 'agTextColumnFilter', editable: true, },
        { field: 'score', filter: 'agNumberColumnFilter', sortable: true, editable: true, },
        { field: 'query_number', editable: true, },
        { field: 'section', filter: 'agTextColumnFilter', editable: true, },
      ],
      summarise: null,
      task_id: null,

      // Filters
      queryText: 'What methods are available to measure anti-mullerian hormone concentrations in young women?',
      queryDate: null,
      queryScore: 0.8,
      queryTypes: new Set(),
    }
  }

  getArticles = (task_id, query_number = 0, interval = 1000) => {
    fetch(variables.API_URL + `/api/ddi_review`,
      {
        headers: {
          'Content-Type': 'application/json;charset=utf-8',
          'Authorization': `Token ${variables.token}`,
        },
      }
    )
      .then(response => {
        console.log(query_number)
        console.log(response.status);
        if (response.ok) {
          return response.json()
        } else {
          ErrorMessage = response.status
          throw Error(response.statusText)
        }
      })
      .then(data => {
        if (data.data === null) {
          this.setState({ loading: true, message: data.message, messageStatus: 202 });
          setTimeout(() => {
            return this.getArticles(task_id, query_number, interval)
          }, interval);
        } else {
          this.setState({
            articles: [...this.state.articles, ...data.data], DetailArticle: data.data[0], loading: false, message: 'Запрос успешно обработан', messageStatus: 200
          });
          if (query_number !== 0) {
            this.state.query_list[query_number - 1].status = 1;
          }
        }
      })
      .catch(error => {
        console.log(error);
        if (ErrorMessage === 500) {
            this.setState({ articles: [], DetailArticle: null, loading: false, message: 'Ошибка сервера', messageStatus: 500 });
        } else if (ErrorMessage === 403) {
            this.setState({ articles: [], DetailArticle: null, loading: false, message: 'Дождитесь окончания предыдушего запроса', messageStatus: 403 });
        } else if (ErrorMessage === 404) {
            this.setState({ articles: [], DetailArticle: null, loading: false, message: 'Сделайте запрос', messageStatus: 202 });
        } else {
            this.setState({ articles: [], DetailArticle: null, loading: false, message: 'Что-то пошло не так', messageStatus: 400 });
        }
        if (query_number !== 0) {
          this.state.query_list[query_number - 1].status = 2;
        }
      })
  }

  createTask() {
    // Отправляем запрос на сервер для получения статей
    this.state.query_list.push({ query: this.state.queryText, status: 0 });
    const query_number = this.state.query_list.length;
    fetch(variables.API_URL + '/api/ddi_review', {
      method: 'POST',
      headers: {
        'Accept': 'application/json',
        'Content-Type': 'application/json;charset=utf-8',
        'Authorization': `Token ${variables.token}`,
      },
      body: JSON.stringify({
        query: this.state.queryText,
        score: this.state.queryScore,
        number_of_query: query_number,
        date: this.state.queryDate,
        type: [...this.state.queryTypes],
      })
    })
      .then(response => {
        console.log(query_number)
        console.log(response.status);
        if (response.ok) {
          return response.json()
        } else {
          ErrorMessage = response.status
          throw Error(response.statusText)
        }
      })
      .then(data => {
        this.setState({
          task_id: data.data,
          message: "Ваш запрос в очереди. Пожайлуста дождитесь результата",
          messageStatus: 201
        });
        this.getArticles(data.data, query_number)
      })
      .catch(error => {
        console.log(error);
        if (ErrorMessage === 500) {
            this.setState({ articles: [], DetailArticle: null, loading: false, message: 'Ошибка сервера', messageStatus: 500 });
        } else if (ErrorMessage === 403) {
            this.setState({ articles: [], DetailArticle: null, loading: false, message: 'Дождитесь окончания предыдушего запроса', messageStatus: 403 });
        } else {
            this.setState({ articles: [], DetailArticle: null, loading: false, message: 'Что-то пошло не так', messageStatus: 400 });
        }
      }
      )
  }

  clearTask() {
    this.setState({ query_list: [], articles: [], DetailArticle: null })
    alert("Таблица очищена!");
  }

  componentDidMount() {
    this.getArticles();
    console.log('start');
  }

  onSelectionChanged = () => {
    const selectedRows = this.gridRef.current.api.getSelectedRows();
    this.setState({ DetailArticle: (selectedRows.length === 1 ? selectedRows[0] : null) })
  }

  changeQueryText = (e) => {
    this.setState({ queryText: e.target.value });
  }

  changeQueryDate = (e) => {
    this.setState({ queryDate: e });
  }

  changeQueryTypes(type) {
    if (this.state.queryTypes.has(type)) {
      this.state.queryTypes.delete(type)
    } else {
      this.state.queryTypes.add(type)
    }
    this.setState({ updateOr: !this.state.updateOr })
  }

  // Summarise

  getSummarise = (task_id, interval = 1000) => {
    fetch(variables.API_URL + `/api/summarise_emb?task_id=${task_id}`,
      {
        headers: {
          'Content-Type': 'application/json;charset=utf-8',
          'Authorization': `Token ${variables.token}`,
        },
      }
    )
      .then((res) => {
        if (res.status == 202) {
          this.setState({ loading: true })
          setTimeout(() => {
            return this.getSummarise(task_id, interval)
          }, interval);
        } else if (res.status == 200) {
          return res.json()
        } else {
          throw Error(res.statusText)
        }
      })
      .then((data) => {
        this.setState({
          summarise: data.data,
          message: 'Суммаризация прошла успешно',
          messageStatus: 200
        });
      })
      .catch((err) => {
        console.log(err);
        this.setState({ message: 'Ошибка при суммаризации', messageStatus: 500, summarise: null })
      });
  }

  createSummariseQuery() {
    fetch(variables.API_URL + '/api/summarise_emb', {
      method: 'POST',
      headers: {
        'Accept': 'application/json',
        'Content-Type': 'application/json;charset=utf-8',
        'Authorization': `Token ${variables.token}`,
      },
      body: JSON.stringify({
        articles: 'some'
      })
    })
      .then((res) => {
        if (res.status == 200) { return res.json() }
        else { throw Error(res.statusText) }
      })
      .then((result) => {
        var task_id = result.data;
        this.setState({ message: 'Отправлено на суммаризацию пожайлуста дождитесь ответа', messageStatus: 201 })
        this.getSummarise(task_id);
      })
      .catch((error) => {
        this.setState({ message: 'Ошибка при суммаризации', messageStatus: 500, summarise: null })
      })
  }

  // MarkUp article

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
        this.setState({ message: 'Произошла ошибка при разметке', loading: false, messageStatus: 500});
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

  suppressCutToClipboard = false;

  onRemoveSelected = () => {

    const selectedData = this.gridRef.current.api.getSelectedRows();
    console.log(selectedData)
    const res = this.gridRef.current.api.applyTransaction({ remove: selectedData });
  }

  onCellValueChanged = (params) => {
    console.log('Callback onCellValueChanged:', params);
    console.log(params.node)
    const res = this.gridRef.current.api.applyTransaction({ remove: [params.node.data] });
  }

  onCutStart = (params) => {
    console.log('Callback onCutStart:', params);
  }

  onCutEnd = (params) => {
    console.log('Callback onCutEnd:', params);
  }

  getRowId = () => {
    return (params) => {
      console.log(params)
      return params.data.code;
    };
  }


  render() {
    const {
      token,
      loading,
      query_list,
      articlesInfo,
      articles,
      DetailArticle,
      message,
      summarise,

      queryText,
      queryDate,
      queryScore,
      allow_page,
      messageStatus,
    } = this.state;

    if (!token) {
      return <Navigate push to="/login" />

    } else if (allow_page === 0) {
      return <Navigate push to="/tematic_review" />
    } else {
      return (
        <>
          <header>
            <nav class="bg-white border-gray-200 px-4 lg:px-6 py-2.5">
              <div class="flex flex-wrap justify-between items-center">
                <div class="flex justify-start items-center">
                  <a href="" class="flex mr-4">
                    <img src="https://flowbite.s3.amazonaws.com/logo.svg" class="mr-3 h-8" alt="FlowBite Logo" />
                    <span class="self-center text-2xl font-semibold whitespace-nowrap">EBM DаtaMed</span>
                  </a>
                  {allow_page === 2?
                      <ul class="flex font-medium flex-row space-x-8">
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
                  <button type="button" class="hidden sm:inline-flex items-center justify-center text-white bg-primary-700 hover:bg-primary-800 focus:ring-4 focus:ring-primary-300 font-medium rounded-lg text-xs px-3 py-1.5 mr-2 dark:bg-primary-600 dark:hover:bg-primary-700 focus:outline-none dark:focus:ring-primary-800"><svg aria-hidden="true" class="mr-1 -ml-1 w-5 h-5" fill="currentColor" viewBox="0 0 20 20" xmlns="http://www.w3.org/2000/svg"><path fill-rule="evenodd" d="M10 5a1 1 0 011 1v3h3a1 1 0 110 2h-3v3a1 1 0 11-2 0v-3H6a1 1 0 110-2h3V6a1 1 0 011-1z" clip-rule="evenodd"></path></svg> Действие</button>
                </div>
              </div>
            </nav>
            <nav class="bg-white border-gray-200 px-6">
              <div class="max-w-screen-xl py-3">
                <div class="flex items-start">
                  <button id="toggleSidebar" aria-expanded="true" aria-controls="sidebar" class="hidden p-2 mr-3 text-gray-600 rounded cursor-pointer lg:inline hover:text-gray-900 hover:bg-gray-100 dark:text-gray-400 dark:hover:text-white dark:hover:bg-gray-700" data-bs-toggle="collapse" data-bs-target="#sidebar" aria-label="Toggle navigation">
                    <svg class="w-6 h-6" fill="currentColor" viewBox="0 0 20 20" xmlns="http://www.w3.org/2000/svg"><path fill-rule="evenodd" d="M3 5a1 1 0 011-1h12a1 1 0 110 2H4a1 1 0 01-1-1zM3 10a1 1 0 011-1h6a1 1 0 110 2H4a1 1 0 01-1-1zM3 15a1 1 0 011-1h12a1 1 0 110 2H4a1 1 0 01-1-1z" clip-rule="evenodd"></path></svg>
                  </button>
                  <label for="topbar-search" class="sr-only">Поисковый запрос</label>
                  <div className='w-full'>
                    <label for="search" class="mb-2 text-sm font-medium text-gray-900 sr-only dark:text-white">Поисковый запрос</label>
                    <div class="relative mt-1 w-full">
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
                        placeholder="Поисковый запрос"
                        value={queryText}
                        onChange={this.changeQueryText}
                        aria-label="Search" />
                      <button type="submit" value="Найти" onClick={() => this.createTask()} class="text-white absolute right-2.5 bottom-2.5 bg-blue-700 hover:bg-blue-800 focus:ring-4 focus:outline-none focus:ring-blue-300 font-medium rounded-lg text-sm px-4 py-2">Найти</button>
                    </div>
                  </div>
                </div>
              </div>
            </nav>
          </header >
          <main>
            <div>
              <div className="container-fluid">
                <div className="row align-items-stretch b-height">
                  <aside id="sidebar" className="h-screen col-md-2 my-3 bg-white collapse show width border rounded-3 g-0">
                    <div className="accordion accordion-flush" id="accordionFlushExample">
                      <div className="accordion-item">
                        <h2 className="accordion-header" id="flush-headingOne">
                          <button className="accordion-button collapsed" type="button" data-bs-toggle="collapse" data-bs-target="#flush-collapseOne" aria-expanded="false" aria-controls="flush-collapseOne">
                            Период поиска
                          </button>
                        </h2>
                        <div id="flush-collapseOne" className="collapse show multi-collapse" aria-labelledby="flush-headingOne" data-bs-target="#accordionFlushExample">
                          <div className="accordion-body">
                            <div class="form-check">
                              <input
                                className="form-check-input"
                                type="radio"
                                name="flexRadioDefault"
                                id="flexRadioDefault1"
                                checked={queryDate === 1}
                                onChange={() => this.changeQueryDate(1)}
                              />
                              <label class="form-check-label" for="flexRadioDefault1">
                                1 год
                              </label>
                            </div>
                            <div class="form-check">
                              <input
                                className="form-check-input"
                                type="radio"
                                name="flexRadioDefault"
                                id="flexRadioDefault2"
                                checked={queryDate === 3}
                                onChange={() => this.changeQueryDate(3)}
                              />
                              <label class="form-check-label" for="flexRadioDefault2">
                                3 года
                              </label>
                            </div>
                            <div class="form-check">
                              <input
                                className="form-check-input"
                                type="radio"
                                name="flexRadioDefault"
                                id="flexRadioDefault3"
                                checked={queryDate === 5}
                                onChange={() => this.changeQueryDate(5)}
                              /><label class="form-check-label" for="flexRadioDefault3">
                                5 лет
                              </label>
                            </div>
                            <div class="form-check">
                              <input
                                className="form-check-input"
                                type="radio"
                                name="flexRadioDefault"
                                id="flexRadioDefault1"
                                checked={!queryDate}
                                onChange={() => this.changeQueryDate(null)}
                              />
                              <label class="form-check-label" for="flexRadioDefault4">
                                            > 5 лет
                              </label>
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
                                checked={this.state.queryTypes.has('Clinical Trial')}
                                onChange={() => this.changeQueryTypes('Clinical Trial')}
                              />
                              <label className="form-check-label" for="CheckboxClinicalTrial">Clinical Trial</label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxMetaAnalysys"
                                name="CheckboxMetaAnalysys"
                                checked={this.state.queryTypes.has('Meta Analysys')}
                                onChange={() => this.changeQueryTypes('Meta Analysys')}
                              />
                              <label className="form-check-label" for="CheckboxMetaAnalysys">Meta Analysys</label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxRandomizedControlledTrial"
                                name="CheckboxRandomizedControlledTrial"
                                checked={this.state.queryTypes.has('Randomized Controlled Trial')}
                                onChange={() => this.changeQueryTypes('Randomized Controlled Trial')}
                              />
                              <label className="form-check-label" for="CheckboxRandomizedControlledTrial">Randomized Controlled Trial</label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxReview"
                                name="CheckboxReview"
                                checked={this.state.queryTypes.has('Review')}
                                onChange={() => this.changeQueryTypes('Review')}
                              />
                              <label className="form-check-label" for="CheckboxReview">Review</label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxSystematicReview"
                                name="CheckboxSystematicReview"
                                checked={this.state.queryTypes.has('Systematic Review')}
                                onChange={() => this.changeQueryTypes('Systematic Review')}
                              />
                              <label className="form-check-label" for="CheckboxSystematicReview">Systematic Review</label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxJournalArticle"
                                name="CheckboxJournalArticle"
                                checked={this.state.queryTypes.has('Journal Article')}
                                onChange={() => this.changeQueryTypes('Journal Article')}
                              />
                              <label className="form-check-label" for="CheckboxJournalArticle">Journal Article</label>
                            </div>

                          </div>
                        </div>
                      </div>
                      <div className="accordion-item">
                        <h2 class="accordion-header" id="flush-headingThree">
                          <button class="accordion-button collapsed" data-target='#flush-collapseThree' type="button" data-bs-toggle="collapse" data-bs-target="#flush-collapseThree" aria-expanded="false" aria-controls="flush-collapseThree">
                            Точность
                          </button>
                        </h2>
                        <div id="flush-collapseThree" className="collapse show multi-collapse" aria-labelledby="flush-headingFour" data-bs-target="#accordionFlushExample">
                          <div className="accordion-body">
                            <p>Требуемая точность = {queryScore.toFixed(2)}</p>
                            <Slider
                              axis="x"
                              x={queryScore}
                              xmax={1}
                              xmin={0}
                              xstep={0.01}
                              onChange={({ x }) => this.setState({ queryScore: x })}
                            />
                          </div>
                        </div>
                      </div>
                      <div className="accordion-item">
                        <h2 class="accordion-header" id="flush-headingThree">
                          <button class="accordion-button collapsed" data-target='#flush-collapseThree' type="button" data-bs-toggle="collapse" data-bs-target="#flush-collapseFive" aria-expanded="false" aria-controls="flush-collapseThree">
                            Выделенные сущности
                          </button>
                        </h2>
                        <div id="flush-collapseFive" className="collapse show multi-collapse" aria-labelledby="flush-headingFour" data-bs-target="#accordionFlushExample">
                          <div className="accordion-body">
                            {Object.entries(obj_color).map(tag => <p class="pb-2 mb-3 border-bottom" style={{ color: `${tag[1]}` }}>{tag[0]}.</p>)}
                          </div>
                        </div>
                      </div>
                    </div>
                  </aside>
                  <section class="col p-3 m-3 border rounded-3 bg-white overflow-auto">
                    <div class="accordion accordion-flush" id="accordion">
                      <div class="accordion-item">
                        <h2 class="accordion-header" id="">
                          {message?
                              messageStatus > 299?
                                <p class="pb-2 mb-3 border-bottom" style={{ color: 'red' }}>{message}.</p>
                                : messageStatus === 200 ?
                                <p class="pb-2 mb-3 border-bottom" style={{ color: 'green' }}>{message}.</p>
                                :
                                <p class="pb-2 mb-3 border-bottom" style={{ color: 'black' }}>{message}.</p>
                              :null
                          }
                          <button class="accordion-button collapsed" type="button" data-bs-toggle="collapse" data-bs-target="#flush-collapseSeven" aria-expanded="false" aria-controls="flush-collapseSeven">
                            Запросы
                          </button>
                        </h2>
                        <div id="flush-collapseSeven" class="collapse multi-collapse" aria-labelledby="flush-headingSeven" data-bs-target="#accordionFlushExample">
                          <div class="accordion-body">
                            {query_list?.map((query, index) =>
                              query.status === 2 ?
                                <p class="pb-2 mb-3 border-bottom" style={{ color: 'red' }}>{index + 1} - {query.query}.</p>
                                : query.status === 1 ?
                                  <p class="pb-2 mb-3 border-bottom" style={{ color: 'green' }}>{index + 1} - {query.query}.</p>
                                  :
                                  <p class="pb-2 mb-3 border-bottom" style={{ color: 'black' }}>{index + 1} - {query.query}.</p>
                            )}
                          </div>
                          <div class="accordion-body">
                            <input className="text-white right-2.5 my-4 bottom-2.5 bg-blue-700 hover:bg-blue-800 focus:ring-4 focus:outline-none focus:ring-blue-300 font-medium rounded-lg text-sm px-4 py-2" type="submit" value="Очистить" onClick={() => this.clearTask()} />
                          </div>
                        </div>
                      </div>
                      <div>
                        {summarise ?
                          <>
                            <p>Summarise</p>
                            <p>{summarise}</p>
                          </>
                          : <input className="text-white right-2.5 my-4 bottom-2.5 bg-blue-700 hover:bg-blue-800 focus:ring-4 focus:outline-none focus:ring-blue-300 font-medium rounded-lg text-sm px-4 py-2" type="submit" value="Суммаризовать" onClick={() => this.createSummariseQuery()} />}
                      </div>
                    </div>
                    <div>
                      <div class="bd-example">
                        <div class="tab-content" id="myTabContent">
                          <div class="tab-pane fade active show" id="home" role="tabpanel" aria-labelledby="home-tab">
                            <div class="container-fluid g-0">
                              <div className="ag-theme-alpine ag-theme-acmecorp" style={{ height: 700 }}>
                                <AgGridReact
                                  ref={this.gridRef}
                                  rowData={articles}
                                  columnDefs={articlesInfo}
                                  pagination={true}
                                  rowSelection={'single'}
                                  onSelectionChanged={this.onSelectionChanged}
                                  suppressCutToClipboard={this.suppressCutToClipboard}
                                  onCellValueChanged={this.onCellValueChanged}
                                  onCutStart={this.onCutStart}

                                  onCutEnd={this.onCutEnd}
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
                    </div>
                  </section>

                  <aside id="sidebar2" class="col-md-4 h-screen collapse show width col p-3 my-3 border rounded-3 bg-white">
                    <h3 class="pb-2 mb-3 border-bottom">Подробное описание</h3>
                    <nav class="small" id="toc">
                      {DetailArticle ?
                        <div class="card mb-3">
                          <div class="card-body">
                            <a href={DetailArticle.url} class="card-title link-primary text-decoration-none h5" target="_blank"> {DetailArticle.titl} </a>
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
                            {summarise ?
                              <>
                                <p>Summarise</p>
                                <p>{summarise}</p>
                              </>
                              : loading ?
                                <p>Loading...</p>
                                : <input className="text-white right-2.5 my-4 bottom-2.5 bg-blue-700 hover:bg-blue-800 focus:ring-4 focus:outline-none focus:ring-blue-300 font-medium rounded-lg text-sm px-4 py-2" type="submit" value="Разметить" onClick={() => this.markUpArticle(DetailArticle)} />
                            }
                            <input className="text-white right-2.5 my-4 bottom-2.5 bg-red-700 hover:bg-red-800 focus:ring-4 focus:outline-none focus:ring-red-300 font-medium rounded-lg text-sm px-4 py-2" type="submit" value="Удалить" onClick={() => this.onRemoveSelected()} />
                          </div>
                        </div>
                        : null}
                    </nav>
                  </aside>

                </div>
              </div>
            </div>
          </main>
        </>
      )
    }
  }
}